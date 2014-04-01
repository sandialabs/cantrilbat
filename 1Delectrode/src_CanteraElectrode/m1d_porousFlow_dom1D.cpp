/**
 * @file m1d_porousFlow_dom1D.cpp
 */

/*
 *   $Id: m1d_porousFlow_dom1D.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */

#include "m1d_porousFlow_dom1D.h"

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

#include <cantera/thermo/StoichSubstance.h>

using namespace std;
using namespace Cantera;

namespace m1d
{
 
  //=====================================================================================================================
  porousFlow_dom1D::porousFlow_dom1D(BulkDomainDescription & bdd) :
    BulkDomain1D(bdd),
    porosity_Cell_(0),
    porosity_Cell_old_(0),
    temp_Curr_(TemperatureReference_)
  {

  }
  //=====================================================================================================================
  porousFlow_dom1D::porousFlow_dom1D(const porousFlow_dom1D &r) :
    BulkDomain1D(r.BDD_),
    porosity_Cell_(0),
    porosity_Cell_old_(0),
   temp_Curr_(TemperatureReference_)
  {
    porousFlow_dom1D::operator=(r);
  }
  //=====================================================================================================================
  porousFlow_dom1D::~porousFlow_dom1D()
  {

  }
  //=====================================================================================================================
  porousFlow_dom1D &
  porousFlow_dom1D::operator=(const porousFlow_dom1D &r)
  {
    if (this == &r) {
      return *this;
    }
    // Call the parent assignment operator
    BulkDomain1D::operator=(r);

    porosity_Cell_ = r.porosity_Cell_;
    porosity_Cell_old_ = r.porosity_Cell_old_;
    temp_Curr_ = r.temp_Curr_;


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
  porousFlow_dom1D::domain_prep(LocalNodeIndices *li_ptr)
  {
    /*
     * First call the parent domain prep to get the node information
     */
    BulkDomain1D::domain_prep(li_ptr);


    double porosity = -1.0;

    porosity_Cell_.resize(NumLcCells, porosity);
    porosity_Cell_old_.resize(NumLcCells, porosity);

  }

  //=====================================================================================================================
} //namespace m1d
//=====================================================================================================================


