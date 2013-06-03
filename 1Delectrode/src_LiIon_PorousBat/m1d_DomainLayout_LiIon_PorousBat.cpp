/**
 * @file m1d_DomainLayout_LiIon_PorousBat.cpp
 *
 */

/*
 *  $Id: m1d_DomainLayout_LiIon_PorousBat.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_DomainLayout_LiIon_PorousBat.h"

#include "m1d_BDT_porSeparator_LiIon.h"
#include "m1d_BDT_porAnode_LiIon.h"
#include "m1d_BDT_porCathode_LiIon.h"
#include "m1d_SDT_AnodeCollector.h"
#include "m1d_SDT_CathodeCollector.h"

#include <iostream>

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace std;

namespace m1d
{

//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat(ProblemStatement *psInput_ptr) :
  DomainLayout(psInput_ptr), ProbNum_(0)
{
  pscInput_ptr_ = dynamic_cast<ProblemStatementCell *>(psInput_ptr);
  if (!pscInput_ptr_) {
    Cantera::CanteraError("DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat()",
			  "Bad dynamic cast");
  }
  InitializeDomainPicture();
}
//===========================================================================
DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat(int probNum, ProblemStatement *psInput_ptr) :
  DomainLayout(), ProbNum_(probNum)
{
  pscInput_ptr_ = dynamic_cast<ProblemStatementCell *>(psInput_ptr);
  if (!pscInput_ptr_) {
    Cantera::CanteraError("DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat()",
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
DomainLayout_LiIon_PorousBat::~DomainLayout_LiIon_PorousBat()
{
}
//====================================================================================================================
DomainLayout_LiIon_PorousBat::DomainLayout_LiIon_PorousBat(const DomainLayout_LiIon_PorousBat &r) :
  DomainLayout(), ProbNum_(0)
{
  *this = r;
}
//====================================================================================================================
DomainLayout_LiIon_PorousBat &
DomainLayout_LiIon_PorousBat::operator=(const DomainLayout_LiIon_PorousBat &r)
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
    if ( !( PSinput.separatorThickness_ > 0.0 ) ) 
      throw CanteraError("DomainLayout_LiIon_PorousBat::malloc_domains()",
			 "separator thickness not specified");
    if ( !( PSinput.anode_input_->electrodeGrossThickness > 0.0 ) ) 
      throw CanteraError("DomainLayout_LiIon_PorousBat::malloc_domains()",
			 "anode thickness not specified");
    if ( !( PSinput.cathode_input_->electrodeGrossThickness > 0.0 ) ) 
      throw CanteraError("DomainLayout_LiIon_PorousBat::malloc_domains()",
			 "cathode thickness not specified");

    double startZ = 0.0;
    double anodeSize = PSinput.anode_input_->electrodeGrossThickness;
    double sepSize = PSinput.separatorThickness_;
    double cathodeSize = PSinput.cathode_input_->electrodeGrossThickness;

    double endZ = startZ + anodeSize;
    BulkDomainDescription *bdd = new BDT_porAnode_LiIon(this);
    // We refine the grid in the anode to get rid of stair step profiles
    //addBulkDomainToRightEnd(bdd, numNodesEach, startZ, endZ);
    int numNodesA = pscInput_ptr_->initDefaultNumCVsAnode_;
    addBulkDomainToRightEnd(bdd, numNodesA, startZ, endZ);

    SDT_AnodeCollector * dirLeft = new SDT_AnodeCollector(this, 1);
    SurfDomainDescription *sddL = dirLeft;
    addSurfDomainToLeftEnd(sddL, bdd);


    SDT_Mixed * anodeSepInterface = new SDT_Mixed(this);
    SurfDomainDescription *sddR = anodeSepInterface;
    addSurfDomainToRightEnd(sddR, bdd);

    startZ = endZ;
    endZ = startZ + sepSize;
    bdd = new BDT_porSeparator_LiIon(this);
    int numNodesS = pscInput_ptr_->initDefaultNumCVsSeparator_;
    addBulkDomainToRightEnd(bdd, numNodesS, startZ, endZ);

    SDT_Mixed * cathodeSepInterface = new SDT_Mixed(this);
    SurfDomainDescription *sddR2 = cathodeSepInterface;
    addSurfDomainToRightEnd(sddR2, bdd);

    startZ = endZ;
    endZ = startZ + cathodeSize;
    bdd = new BDT_porCathode_LiIon(this);
    int numNodesC = pscInput_ptr_->initDefaultNumCVsCathode_;
    addBulkDomainToRightEnd(bdd, numNodesC, startZ, endZ);

    SDT_CathodeCollector * cc = new SDT_CathodeCollector(this, 1);
    SurfDomainDescription *sddR3 = cc;
    addSurfDomainToRightEnd(sddR3, bdd);
  }
//===========================================================================
}
