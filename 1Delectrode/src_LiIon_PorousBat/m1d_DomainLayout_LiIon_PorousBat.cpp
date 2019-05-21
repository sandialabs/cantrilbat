/**
 * @file m1d_DomainLayout_LiIon_PorousBat.cpp
 *
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "LiIon_PorousBat.h"

#include "m1d_DomainLayout_LiIon_PorousBat.h"
#include "m1d_BDD_porSeparator_LiIon.h"
#include "m1d_BDD_porAnode_LiIon.h"
#include "m1d_BDD_porCathode_LiIon.h"
#include "m1d_SDD_AnodeCollector.h"
#include "m1d_SDD_CathodeCollector.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_BDD_porousElectrode.h"
#include "m1d_SDD_ElectrodeSepInterface.h"

#include <iostream>

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

namespace m1d
{
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat(ProblemStatement* psInput_ptr) :
    DomainLayout(psInput_ptr), 
    ProbNum_(0)
{
    pscInput_ptr_ = dynamic_cast<ProblemStatementCell*>(psInput_ptr);
    if (!pscInput_ptr_) {
        Zuzax::ZuzaxError("DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat()",
                              "Bad dynamic cast");
    }
}
//===========================================================================
DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat(int probNum, ProblemStatement* psInput_ptr) :
    DomainLayout(), ProbNum_(probNum)
{
    pscInput_ptr_ = dynamic_cast<ProblemStatementCell*>(psInput_ptr);
    if (!pscInput_ptr_) {
        Zuzax::ZuzaxError("DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat()",
                              "Bad dynamic cast");
    }
    if (probNum == 1) {
        /*
         *  Layout the problem. This call will call malloc_domains() below.
         */
        InitializeDomainPicture();
    } else {
        std::cerr << "DomainLayout constructor: Unknown problem # " << probNum << endl;
        exit(-1);
    }
}
//====================================================================================================================
DomainLayout_LiIon_PorousBat::~DomainLayout_LiIon_PorousBat()
{
}
//====================================================================================================================
DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat(const DomainLayout_LiIon_PorousBat& r) :
    DomainLayout(), ProbNum_(0)
{
    *this = r;
}
//====================================================================================================================
DomainLayout_LiIon_PorousBat&
DomainLayout_LiIon_PorousBat::operator=(const DomainLayout_LiIon_PorousBat& r)
{
    if (this == &r) {
        return *this;
    }
    DomainLayout::operator=(r);

    ProbNum_ = r.ProbNum_;

    return *this;
}
//====================================================================================================================
// Allocate the domain structure
void
DomainLayout_LiIon_PorousBat::malloc_domains()
{

    /*
     *  Here we lay out the size of the domain
     */
    // First check to see if we have thickness for each layer
    if (!(PSinput.separatorThickness_ > 0.0))
        throw ZuzaxError("DomainLayout_LiIon_PorousBat::malloc_domains()",
                           "separator thickness not specified");
   /*
    if (!(PSinput.anode_input_->electrodeGrossThickness > 0.0)) 
        throw ZuzaxError("DomainLayout_LiIon_PorousBat::malloc_domains()",
                           "anode thickness not specified");
    if (!(PSinput.cathode_input_->electrodeGrossThickness > 0.0))
        throw ZuzaxError("DomainLayout_LiIon_PorousBat::malloc_domains()",
                           "cathode thickness not specified");
   */
    double startZ = 0.0;
    double anodeSize = PSinput.anode_input_->electrodeGrossThickness;
    double sepSize = PSinput.separatorThickness_;
    double cathodeSize = PSinput.cathode_input_->electrodeGrossThickness;

    double endZ = startZ + anodeSize;


    //ProblemStatementCell* psc_ptr = &PSinput;
    //ELECTRODE_KEY_INPUT* ai = psc_ptr->anode_input_;
    //ELECTRODE_KEY_INPUT* ci = psc_ptr->cathode_input_;


    BulkDomainDescription* bdd = new BDD_porAnode_LiIon(this, "PorousLiIonAnode");
    // We refine the grid in the anode to get rid of stair step profiles
    //addBulkDomainToRightEnd(bdd, numNodesEach, startZ, endZ);
    int numNodesA = pscInput_ptr_->initDefaultNumCVsAnode_;
    addBulkDomainToRightEnd(bdd, numNodesA, startZ, endZ);

    SDD_AnodeCollector* dirLeft = new SDD_AnodeCollector(this, 1);
    SurfDomainDescription* sddL = dirLeft;
    addSurfDomainToLeftEnd(sddL, bdd);

    SDD_ElectrodeSepInterface* anodeSepInterface = new SDD_ElectrodeSepInterface(this);
    SurfDomainDescription* sddR = anodeSepInterface;
    addSurfDomainToRightEnd(sddR, bdd);

    startZ = endZ;
    endZ = startZ + sepSize;
    bdd = new BDD_porSeparator_LiIon(this);
    int numNodesS = pscInput_ptr_->initDefaultNumCVsSeparator_;
    addBulkDomainToRightEnd(bdd, numNodesS, startZ, endZ);

    SDD_ElectrodeSepInterface* cathodeSepInterface = new SDD_ElectrodeSepInterface(this);
    SurfDomainDescription* sddR2 = cathodeSepInterface;
    addSurfDomainToRightEnd(sddR2, bdd);

    startZ = endZ;
    endZ = startZ + cathodeSize;
    bdd = new BDD_porCathode_LiIon(this);
    int numNodesC = pscInput_ptr_->initDefaultNumCVsCathode_;
    addBulkDomainToRightEnd(bdd, numNodesC, startZ, endZ);

    SDD_CathodeCollector* cc = new SDD_CathodeCollector(this, 1);
    SurfDomainDescription* sddR3 = cc;
    addSurfDomainToRightEnd(sddR3, bdd);
}
//=====================================================================================================
}
