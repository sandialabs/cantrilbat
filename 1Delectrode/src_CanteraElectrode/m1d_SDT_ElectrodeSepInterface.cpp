/**
 * @file SDT_ElectrodeSepInterface.cpp
 * Implementation of an interface between anode-separator or cathode-separator
 * (see class \link m1d::SDT_ElectrodeSepInterface SDT_ElectrodeSepInterface\endlink).
 */

/*
 * Copywrite 2014 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_SDT_ElectrodeSepInterface.h"
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
SDT_ElectrodeSepInterface::SDT_ElectrodeSepInterface(DomainLayout *dl_ptr, const char *domainName) :
    SDT_Mixed(dl_ptr, domainName)
{
}
//==================================================================================================================================
SDT_ElectrodeSepInterface::SDT_ElectrodeSepInterface(const SDT_ElectrodeSepInterface &r) :
    SDT_Mixed(r.DL_ptr_, r.DomainName)
{
  *this = r;
}
//==================================================================================================================================
SDT_ElectrodeSepInterface::~SDT_ElectrodeSepInterface()
{
}
//==================================================================================================================================
SDT_ElectrodeSepInterface &
SDT_ElectrodeSepInterface::operator=(const SDT_ElectrodeSepInterface &r)
{
  if (this == &r) {
    return *this;
  }

  SDT_Mixed::operator=(r);

  return *this;
}
//==================================================================================================================================
} /* End of Namespace */
//==================================================================================================================================
