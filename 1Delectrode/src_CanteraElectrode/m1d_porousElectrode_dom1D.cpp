/**
 * @file m1d_porousElectrode_dom1D.cpp
 */

/*
 *   $Id: m1d_porousElectrode_dom1D.cpp 564 2013-03-08 23:35:51Z hkmoffa $
 */

#include "m1d_porousElectrode_dom1D.h"

using namespace std;
using namespace Cantera;

namespace m1d
{
 
  //=====================================================================================================================
  porousElectrode_dom1D::porousElectrode_dom1D(BulkDomainDescription & bdd) :
      porousFlow_dom1D(bdd),
      maxElectrodeSubIntegrationSteps_(0)
  {

  }
  //=====================================================================================================================
  porousElectrode_dom1D::porousElectrode_dom1D(const porousElectrode_dom1D &r) :
    porousFlow_dom1D(r.BDD_)
  {
    porousElectrode_dom1D::operator=(r);
  }
  //=====================================================================================================================
  porousElectrode_dom1D::~porousElectrode_dom1D()
  {

  }
  //=====================================================================================================================
  porousElectrode_dom1D &
  porousElectrode_dom1D::operator=(const porousElectrode_dom1D &r)
  {
    if (this == &r) {
      return *this;
    }
    // Call the parent assignment operator
    porousFlow_dom1D::operator=(r);

    maxElectrodeSubIntegrationSteps_ = r.maxElectrodeSubIntegrationSteps_;

    return *this;
  }
  //=====================================================================================================================
  // Prepare all of the indices for fast calculation of the residual
  /*
   *  Ok, at this point, we will have figured out the number of equations
   *  to be calculated at each node point. The object NodalVars will have
   *  been fully formed.
   *
   *  We use this to figure out what local node numbers/ cell numbers are
   *  needed and to set up indices for their efficient calling.
   *
   *  Child objects of this one will normally call this routine in a
   *  recursive fashion.
   */
  void
  porousElectrode_dom1D::domain_prep(LocalNodeIndices *li_ptr)
  {
    /*
     * First call the parent domain prep to get the node information
     */
    porousFlow_dom1D::domain_prep(li_ptr);

  }
  //====================================================================================================================
  //  An electrode object must be created and initialized for every cell in the domain
  /*
   * Create electrode objects for every cell.  
   * Correct the volume and number of moles of 
   * active material within each of these electrode 
   * objects to correspond to the discretized volume.
   */
  void
  porousElectrode_dom1D::instantiateElectrodeCells() 
  {
  }
//====================================================================================================================
int porousElectrode_dom1D::getMaxSubGridTimeSteps() const
{
    return maxElectrodeSubIntegrationSteps_;
}
  //=====================================================================================================================
} //namespace m1d
//=====================================================================================================================

