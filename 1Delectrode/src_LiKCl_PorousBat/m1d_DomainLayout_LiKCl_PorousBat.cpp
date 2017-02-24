/**
 * @file m1d_DomainLayout_LiKCl_infPorousBat.cpp
 *
 */

#include "m1d_defs.h"
#include "m1d_BulkDomainTypes.h"
#include "m1d_SDD_Mixed.h"

#include "m1d_DomainLayout_LiKCl_PorousBat.h"
#include "m1d_BDT_porousLiKCl.h"
#include "m1d_BDT_porAnode_LiKCl.h"
#include "m1d_BDT_porCathode_LiKCl.h"
#include "m1d_SDD_AnodeCollector.h"
#include "m1d_SDD_CathodeCollector.h"

#include "m1d_exception.h"

#include <iostream>

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

namespace m1d
{
//=====================================================================================================================
DomainLayout_LiKCl_infPorousBat::DomainLayout_LiKCl_infPorousBat(ProblemStatement *psInput_ptr) :
  DomainLayout(psInput_ptr), ProbNum_(0)
{
  pscInput_ptr_ = dynamic_cast<ProblemStatementCell *>(psInput_ptr);
  if (!pscInput_ptr_) {
    ZZCantera::CanteraError("DomainLayout_LiKCl_infPorousBat::DomainLayout_LiKCl_infPorousBat()",
			  "Bad dynamic cast");
  }
  InitializeDomainPicture();
}
//===========================================================================
DomainLayout_LiKCl_infPorousBat::DomainLayout_LiKCl_infPorousBat(int probNum, ProblemStatement *psInput_ptr) :
  DomainLayout(), ProbNum_(probNum)
{
  pscInput_ptr_ = dynamic_cast<ProblemStatementCell *>(psInput_ptr);
  if (!pscInput_ptr_) {
    ZZCantera::CanteraError("DomainLayout_LiKCl_infPorousBat::DomainLayout_LiKCl_infPorousBat()",
			  "Bad dynamic cast");
  }
  if (probNum == 1) {
    /*
     *  Layout the problem. This call will call malloc_domains() below.
     */
    InitializeDomainPicture();
  } else {
    cerr << "DomainLayout constructor: Unknown problem # " << probNum << endl;
    exit(-1);
  }
}
//====================================================================================================================
DomainLayout_LiKCl_infPorousBat::~DomainLayout_LiKCl_infPorousBat()
{
}
//====================================================================================================================
DomainLayout_LiKCl_infPorousBat::DomainLayout_LiKCl_infPorousBat(const DomainLayout_LiKCl_infPorousBat &r) :
  DomainLayout(), ProbNum_(0)
{
  *this = r;
}
//====================================================================================================================
DomainLayout_LiKCl_infPorousBat &
DomainLayout_LiKCl_infPorousBat::operator=(const DomainLayout_LiKCl_infPorousBat &r)
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
  DomainLayout_LiKCl_infPorousBat::malloc_domains()
  {

    /*
     *  Here we lay out the size of the domain
     */
    // First check to see if we have thickness for each layer
    if ( !( PSinput.separatorThickness_ > 0.0 ) ) 
      throw CanteraError("DomainLayout_LiKCl_infPorousBat::malloc_domains()",
			 "separator thickness not specified");
    if (!( PSinput.anode_input_->electrodeGrossThickness > 0.0 ) ) 
      throw CanteraError("DomainLayout_LiKCl_infPorousBat::malloc_domains()",
			 "anode thickness not specified");
    if ( !( PSinput.cathode_input_->electrodeGrossThickness > 0.0 ) ) 
      throw CanteraError("DomainLayout_LiKCl_infPorousBat::malloc_domains()",
			 "cathode thickness not specified");

    double startZ = 0.0;
    double anodeSize = PSinput.anode_input_->electrodeGrossThickness;
    double sepSize = PSinput.separatorThickness_;
    double cathodeSize = PSinput.cathode_input_->electrodeGrossThickness;

    double endZ = startZ + anodeSize;
    BulkDomainDescription *bdd = new BDT_porAnode_LiKCl(this);
    // We refine the grid in the anode to get rid of stair step profiles
    //addBulkDomainToRightEnd(bdd, numNodesEach, startZ, endZ);
    int numNodesA = pscInput_ptr_->initDefaultNumCVsAnode_;
    addBulkDomainToRightEnd(bdd, numNodesA, startZ, endZ);

    SDD_AnodeCollector * dirLeft = new SDD_AnodeCollector(this, 1);
    SurfDomainDescription *sddL = dirLeft;
    addSurfDomainToLeftEnd(sddL, bdd);


    SDD_Mixed * anodeSepInterface = new SDD_Mixed(this);
    SurfDomainDescription *sddR = anodeSepInterface;
    addSurfDomainToRightEnd(sddR, bdd);

    startZ = endZ;
    endZ = startZ + sepSize;
    bdd = new BDT_porousLiKCl(this);
    int numNodesS = pscInput_ptr_->initDefaultNumCVsSeparator_;
    addBulkDomainToRightEnd(bdd, numNodesS, startZ, endZ);

    SDD_Mixed * cathodeSepInterface = new SDD_Mixed(this);
    SurfDomainDescription *sddR2 = cathodeSepInterface;
    addSurfDomainToRightEnd(sddR2, bdd);

    startZ = endZ;
    endZ = startZ + cathodeSize;
    bdd = new BDT_porCathode_LiKCl(this);
    int numNodesC = pscInput_ptr_->initDefaultNumCVsCathode_;
    addBulkDomainToRightEnd(bdd, numNodesC, startZ, endZ);

    SDD_CathodeCollector * cc = new SDD_CathodeCollector(this, 1);
    SurfDomainDescription *sddR3 = cc;
    addSurfDomainToRightEnd(sddR3, bdd);
  }
//===========================================================================
}
