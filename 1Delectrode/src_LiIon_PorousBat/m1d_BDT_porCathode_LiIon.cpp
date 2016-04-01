/**
 * @file m1d_BDT_porCathode_LiIon.cpp
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.x
 */

#include "m1d_BDT_porCathode_LiIon.h"
#include "m1d_porousLiIon_Cathode_dom1D.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_exception.h"
#include "m1d_defs.h"

#include "Electrode.h"
#include "Electrode_Factory.h"

#include "cantera/thermo.h"


extern m1d::ProblemStatementCell PSinput;

using namespace std;

namespace m1d
{

//=====================================================================================================================
BDT_porCathode_LiIon::BDT_porCathode_LiIon(DomainLayout *dl_ptr, std::string domainName) :
    BDD_porousElectrode(dl_ptr, 1, "cathode", domainName),
    m_position(1),
    Capacity_initial_(0.0),
    CapacityDischarged_initial_(0.0),
    CapacityLeft_initial_(0.0)
{

    IsAlgebraic_NE.resize(7,0);
    IsArithmeticScaled_NE.resize(7,0);   
}
//=====================================================================================================================
BDT_porCathode_LiIon::BDT_porCathode_LiIon(const BDT_porCathode_LiIon& r) :
    BDD_porousElectrode(r),
    m_position(1),
    Capacity_initial_(0.0),
    CapacityDischarged_initial_(0.0),
    CapacityLeft_initial_(0.0)
{
    *this = r;
}
//=====================================================================================================================
BDT_porCathode_LiIon::~BDT_porCathode_LiIon()
{ 
}
//=====================================================================================================================
BDT_porCathode_LiIon&
BDT_porCathode_LiIon::operator=(const BDT_porCathode_LiIon& r)
{
    if (this == &r) {
        return *this;
    }

    BDD_porousElectrode::operator=(r);

    m_position = r.m_position;
    Capacity_initial_ = r.Capacity_initial_;
    CapacityDischarged_initial_ = r.CapacityDischarged_initial_;
    CapacityLeft_initial_ = r.CapacityLeft_initial_;

    return *this;
}
//=====================================================================================================================
void
BDT_porCathode_LiIon::ReadModelDescriptions()
{
    BDD_porousElectrode::ReadModelDescriptions();

    Capacity_initial_ = Electrode_->capacity();
    CapacityDischarged_initial_ = Electrode_->capacityDischarged();
    CapacityLeft_initial_ = Electrode_->capacityLeft();
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D*
BDT_porCathode_LiIon::mallocDomain1D()
{
    BulkDomainPtr_ = new porousLiIon_Cathode_dom1D(this);
    return BulkDomainPtr_;
}
//=====================================================================================================================
void
BDT_porCathode_LiIon::DetermineConstitutiveModels()
{
    BDD_porousElectrode::DetermineConstitutiveModels();
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
