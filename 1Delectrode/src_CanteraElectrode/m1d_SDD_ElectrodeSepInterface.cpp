/**
 * @file SDD_ElectrodeSepInterface.cpp
 * Implementation of an interface between anode-separator or cathode-separator
 * (see class \link m1d::SDD_ElectrodeSepInterface SDD_ElectrodeSepInterface\endlink).
 */

/*
 * Copywrite 2014 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_SDD_ElectrodeSepInterface.h"
#include "m1d_SurDomain_CathodeCollector.h"
#include "m1d_exception.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_BC_Battery.h"
#include "m1d_CanteraElectrodeGlobals.h"

#include "m1d_BDD_porousElectrode.h"

//==================================================================================================================================
namespace m1d
{
//==================================================================================================================================
SDD_ElectrodeSepInterface::SDD_ElectrodeSepInterface(DomainLayout *dl_ptr, const char *domainName) :
    SDD_Mixed(dl_ptr, domainName)
{
}
//==================================================================================================================================
SDD_ElectrodeSepInterface::SDD_ElectrodeSepInterface(const SDD_ElectrodeSepInterface &r) :
    SDD_Mixed(r.DL_ptr_, r.DomainName_)
{
    operator=(r);
}
//==================================================================================================================================
SDD_ElectrodeSepInterface::~SDD_ElectrodeSepInterface()
{
}
//==================================================================================================================================
SDD_ElectrodeSepInterface &
SDD_ElectrodeSepInterface::operator=(const SDD_ElectrodeSepInterface &r)
{
  if (this == &r) {
    return *this;
  }

  SDD_Mixed::operator=(r);

  return *this;
}
//==================================================================================================================================
} /* End of Namespace */
//==================================================================================================================================
