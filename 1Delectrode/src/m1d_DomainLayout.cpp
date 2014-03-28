/**
 * @file m1d_DomainLayout.cpp
 *
 */

/*
 *  $Id: m1d_DomainLayout.cpp 540 2013-02-27 22:18:26Z hkmoffa $
 */
#include "m1d_defs.h"
#include "m1d_DomainLayout.h"
#include "m1d_BulkDomainTypes.h"
#include "m1d_SurfDomainTypes.h"
#include "m1d_ProblemStatement.h"

#include "m1d_exception.h"
#include "m1d_globals.h"

#include <iostream>

using namespace std;

namespace m1d
{
//===========================================================================
// Pointer to the global DomainLayout object
//DomainLayout *DL_Global_ptr = 0;
//===========================================================================
DomainLayout::DomainLayout(ProblemStatement *psInput_ptr) :
  NumDomains(0), NumBulkDomains(0), NumSurfDomains(0), NumGbNodes(0), 
  StartXLoc_Domain(0), EndXLoc_Domain(0),
  XLoc_LeftBoundary(0.0), XLoc_RightBoundary(0.0),
  problemResid_(0), psInput_ptr_(psInput_ptr),
  domainList_(0)
{
   m1d::readEnvironmentalVariables();
}
//===========================================================================
  DomainLayout::DomainLayout(std::vector<std::string>  domainList, ProblemStatement *psInput_ptr,
		std::map<std::string,Cantera::MultiPhase*> bulkMap,
		std::map<std::string,int> surfaceMap ) :
  NumDomains(0), NumBulkDomains(0), NumSurfDomains(0), NumGbNodes(0), 
  StartXLoc_Domain(0), EndXLoc_Domain(0),
  XLoc_LeftBoundary(0.0), XLoc_RightBoundary(0.0), 
  problemResid_(0), psInput_ptr_(psInput_ptr),
  domainList_(domainList)
{
}
//===========================================================================
DomainLayout::~DomainLayout() {
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    safeDelete(BulkDomainDesc_global[ibd]);
    safeDelete(BulkDomain1D_List[ibd]);
  }

  for (int isd = 0; isd < NumSurfDomains; isd++) {
    safeDelete(SurfDomainDesc_global[isd]);
    safeDelete(SurDomain1D_List[isd]);
  }
}
//===========================================================================
DomainLayout::DomainLayout(const DomainLayout &r) :
  NumDomains(0), NumBulkDomains(0), NumSurfDomains(0), NumGbNodes(0), 
  StartXLoc_Domain(0), EndXLoc_Domain(0),
  XLoc_LeftBoundary(0.0), XLoc_RightBoundary(0.0),
  problemResid_(0), psInput_ptr_(0),
  domainList_(r.domainList_)
{
  *this = r;
}
//===========================================================================
DomainLayout &
DomainLayout::operator=(const DomainLayout &r)
{
  if (this == &r)
    return *this;

  NumDomains = r.NumDomains;
  NumBulkDomains = r.NumBulkDomains;
  NumSurfDomains = r.NumSurfDomains;
  NumGbNodes = r.NumGbNodes;

  BulkDomainDesc_global = r.BulkDomainDesc_global;
  BulkDomain1D_List = r.BulkDomain1D_List;

  SurfDomainDesc_global = r.SurfDomainDesc_global;
  SurDomain1D_List = r.SurDomain1D_List;

  StartGBNode_Domain = r.StartGBNode_Domain;
  EndGBNode_Domain = r.EndGBNode_Domain;

  locGBNode_SurfDomain_List = r.locGBNode_SurfDomain_List;
  XLoc_LeftBoundary = r.XLoc_LeftBoundary;
  XLoc_RightBoundary = r.XLoc_RightBoundary;

  problemResid_ = r.problemResid_;
  psInput_ptr_  = r.psInput_ptr_;
  domainList_ = r.domainList_;

  printf("not implemented\n");
  exit(-1);

  return *this;
}
//===========================================================================
// Overarching structure for getting DomainLayout and
// DomainDescription Objects initialized and filled.
/*
 *
 */
void
DomainLayout::InitializeDomainPicture()
{

  /*
   *   Get the initial DomainDescription objects malloced
   *    We will invoke calls to addBulkDomain() and addSurfDomain()
   */
  malloc_domains();
  /*
   *  Discover the ordering of the domains and tell the domains what's
   *  next to them.
   */
  updateAdjInfoInDomainDescriptions();
  /*
   *  Tell the domain descriptions what their global nodes are
   */
  updateGbNodeInfoInDomainDescriptions();
  /*
   *  Tell the domain descriptions what their domain boundary positions are
   */
  updateXposInDomainDescriptions();

  /*
   *  Discover what equations are defined on which domains.
   *  - determine how tie conditions are handled between bulk domains
   *    and surface domains.
   */
  SetEqnDescriptions();
}
//===========================================================================
// Allocate the domain structure
void
DomainLayout::malloc_domains()
{
  /*
   *  Calls to addBulkDomainToRightEnd() to set up the initial bulk domains
   *    -> this is
   */

}
//======================================================================================================================
void
DomainLayout::addBulkDomainToRightEnd(BulkDomainDescription *bdd, int numNodes, double startZ, double endZ)
{
  NumDomains++;
  NumBulkDomains++;

  /*
   * Resize the bulk domains lists within the object
   */
  BulkDomainDesc_global.resize(NumBulkDomains);
  BulkDomain1D_List.resize(NumBulkDomains);
  StartGBNode_Domain.resize(NumBulkDomains);
  EndGBNode_Domain.resize(NumBulkDomains);
  StartXLoc_Domain.resize(NumBulkDomains);
  EndXLoc_Domain.resize(NumBulkDomains);

  int iDom = NumBulkDomains - 1;
  /*
   *  Add the domain description to the list
   */
  BulkDomainDesc_global[iDom] = bdd;
  bdd->setID(iDom);
  /*
   * Figure out the starting and ending global node number of the new domain
   */
  if (NumBulkDomains == 1) {
    StartGBNode_Domain[0] = 0;
  } else {
    StartGBNode_Domain[iDom] = EndGBNode_Domain[iDom - 1];
  }
  EndGBNode_Domain[iDom] = StartGBNode_Domain[iDom] + numNodes - 1;
  NumGbNodes = EndGBNode_Domain[iDom] + 1;
  /*
   * Figure out the starting and ending locations of each of the domains
   */
  if (NumBulkDomains == 1) {
    StartXLoc_Domain[iDom] = startZ;
    EndXLoc_Domain[iDom] = endZ;
    XLoc_LeftBoundary = StartXLoc_Domain[iDom];
    XLoc_RightBoundary = EndXLoc_Domain[iDom];
  } else {
    double endLast = EndXLoc_Domain[iDom - 1];
    /*
     * We throw an error if the domains don't stack up against each other
     * at the moment until we figure out a situation where their
     * not stacking up makes sense.
     */
    if (fabs(startZ - endLast) > 1.0E-6) {
      throw m1d_Error("DomainLayout::addBulkDomainToRightEnd", "Node connectivity issue");
    }
    StartXLoc_Domain[iDom] = endLast;
    EndXLoc_Domain[iDom] = endZ;
    XLoc_RightBoundary = endZ;
  }
  /*
   * Now that we have figured out the boundaries, propagate that down to the bulk domain
   * description object.
   */
  bdd->setXposBounds(StartXLoc_Domain[iDom], EndXLoc_Domain[iDom]);
  /*
   *  Attach the bulk domain to the surface domain to the left
   */
  if (NumSurfDomains >= 2) {
    SurfDomainDescription *sddLast = SurfDomainDesc_global[NumSurfDomains-1];
    if (sddLast->RightBulk == 0) {
      sddLast->RightBulk = bdd;
      sddLast->RightDomain = bdd;
      bdd->LeftSurf = sddLast;
    }
  }
}
//===========================================================================
//  Tell the domain descriptions what they need to know
void
DomainLayout::InfoToDomainDescriptions()
{
  updateAdjInfoInDomainDescriptions();

  updateGbNodeInfoInDomainDescriptions();
}
//===========================================================================
//  Tell the domain descriptions what they need to know
/*
 *   Note, that this routine is now basically an internal check routine.
 *   The functionality has been placed elsewhere
 */
void
DomainLayout::updateAdjInfoInDomainDescriptions()
{
  /*
   * Tell each bulk domain description what surface domain is to the
   * left and right of it.
   * We are relying here on a simple alternating layout. This may
   * be relaxed in the future.
   */
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    SurfDomainDescription *leftSurf = SurfDomainDesc_global[ibd];
    SurfDomainDescription *rightSurf = SurfDomainDesc_global[ibd + 1];
    if (leftSurf != bdd->LeftSurf) {
      throw m1d_Error("DomainLayout::updateAdjInfoInDomainDescriptions()",
		      "surface domains don't match as they should");
    }
    if (rightSurf != bdd->RightSurf) {
      throw m1d_Error("DomainLayout::updateAdjInfoInDomainDescriptions()",
		      "surface domains don't match as they should");
    }
    // bdd->setAdjSurfDomains(leftSurf, rightSurf);
  }

  /*
   * Tell each surface domain description what bulk domain is to the
   * left and right of it.
   * We are relying here on a simple alternating layout. This may
   * be relaxed in the future.
   */
  for (int isd = 0; isd < NumSurfDomains; isd++) {
    SurfDomainDescription *sdd = SurfDomainDesc_global[isd];
    sdd->setID(isd);
    BulkDomainDescription *leftBulk = 0;
    if (isd - 1 >= 0) {
      leftBulk = BulkDomainDesc_global[isd - 1];
    }
    BulkDomainDescription *rightBulk = 0;
    if (isd < NumBulkDomains) {
      rightBulk = BulkDomainDesc_global[isd];
    }
    if (leftBulk != sdd->LeftBulk) {
      throw m1d_Error("DomainLayout::updateAdjInfoInDomainDescriptions()",
		      "bulk domains don't match as they should");
    }
    if (rightBulk != sdd->RightBulk) {
      throw m1d_Error("DomainLayout::updateAdjInfoInDomainDescriptions()",
		      "surface domains don't match as they should");
    }
    //sdd->setAdjBulkDomains(leftBulk, rightBulk);
  }

}

//===========================================================================
//  Tell the domain descriptions what they need to know
void
DomainLayout::updateGbNodeInfoInDomainDescriptions()
{
  /*
   * Tell each bulk domain description what global nodes it encompasses
   */
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    int leftGbNode = StartGBNode_Domain[ibd];
    int rightGbNode = EndGBNode_Domain[ibd];
    bdd->setGbNodeBounds(leftGbNode, rightGbNode);
  }

  for (int isd = 0; isd < NumSurfDomains; isd++) {
    SurfDomainDescription *sdd = SurfDomainDesc_global[isd];
    int gbn = locGBNode_SurfDomain_List[isd];
    sdd->setGbNode(gbn);
  }
}
//===========================================================================
// Add a surface domain to the left end
/*
 *  @param  sddL  Left Surface domain description
 *   @param bdd    bulk domain description
 */
void
DomainLayout::addSurfDomainToLeftEnd(SurfDomainDescription *sddL, BulkDomainDescription *bddA)
{
  int fbd = -1;
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    if (bdd == bddA) {
      fbd = ibd;
    }
  }
  if (fbd == -1) {
    throw m1d_Error("", "");
  }
  int id = NumSurfDomains;
  NumDomains++;
  NumSurfDomains++;
  SurfDomainDesc_global.push_back(sddL);
  bddA->LeftSurf = sddL;
  sddL->setID(id);

  BulkDomainDescription *bddL = 0;
  if (fbd > 0) {
    bddL = BulkDomainDesc_global[fbd - 1];
    bddL->RightSurf = sddL;
  }
  sddL->setAdjBulkDomains(bddL, bddA);
  SurDomain1D_List.resize(NumSurfDomains);
  locGBNode_SurfDomain_List.push_back(StartGBNode_Domain[fbd]);
}
//======================================================================================================================
// Add a surface domain to the right end of the a bulk domain.
/*
 *  @param  sddr    Surface domain description to be added to the right end
 *   @param bddL    bulk domain description that is to the left of the current surface domain
 */
void
DomainLayout::addSurfDomainToRightEnd(SurfDomainDescription *sddR, BulkDomainDescription *bddL)
{
  int fbd = -1;
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    if (bdd == bddL) {
      fbd = ibd;
    }
  }
  if (fbd == -1) {
    throw m1d_Error("DomainLayout::addSurfDomainToRightEnd", 
		    "Couldn't find left bulk domain in the list of bulk domains.");
  }
  int id = NumSurfDomains;
  NumDomains++;
  NumSurfDomains++;
  SurfDomainDesc_global.push_back(sddR);
  bddL->RightSurf = sddR;
  sddR->setID(id);

  BulkDomainDescription *bddR = 0;
  if (fbd < NumBulkDomains - 1) {
    bddR = BulkDomainDesc_global[fbd + 1];
    bddR->LeftSurf = sddR;
  }
  sddR->setAdjBulkDomains(bddL, bddR);

  SurDomain1D_List.resize(NumSurfDomains);
  locGBNode_SurfDomain_List.push_back(EndGBNode_Domain[fbd]);
}
//======================================================================================================================
// Tell the domain descriptions what their domain boundary positions are
void
DomainLayout::updateXposInDomainDescriptions()
{
  double xleft, xright;
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    xleft = StartXLoc_Domain[ibd];
    xright = EndXLoc_Domain[ibd];
    bdd->setXposBounds(xleft, xright);
  }
}
//======================================================================================================================
void
DomainLayout::InitializeXposNodes(GlobalIndices *gi_ptr)
{
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    bdd->InitializeXposNodes(gi_ptr);
  }
}
//======================================================================================================================
void
DomainLayout::SetEqnDescriptions()
{
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    bdd->SetEquationDescription();
  }

  for (int isd = 0; isd < NumSurfDomains; isd++) {
    SurfDomainDescription *sdd = SurfDomainDesc_global[isd];
    sdd->SetEquationDescription();
  }
}
//======================================================================================================================
void
DomainLayout::generateDomain1D(ProblemResidEval *problemResid_ptr)
{
  problemResid_ = problemResid_ptr;

  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomainDescription *bdd = BulkDomainDesc_global[ibd];
    BulkDomain1D_List[ibd] = bdd->mallocDomain1D();
  }

  for (int isd = 0; isd < NumSurfDomains; isd++) {
    SurfDomainDescription *sdd = SurfDomainDesc_global[isd];
    SurDomain1D_List[isd] = sdd->mallocDomain1D();
  }
}
//======================================================================================================================
void
DomainLayout::set_TP_Reference(double temp_ref, double pres_ref)
{
  if (!BulkDomain1D_List[0]) {
    throw m1d_Error("DomainLayout::setTemperatureReference",
                    "Need to malloc Domain1D structures first");
  }
  for (int ibd = 0; ibd < NumBulkDomains; ibd++) {
    BulkDomain1D_List[ibd]->TemperatureReference_ = temp_ref;
    BulkDomain1D_List[ibd]->PressureReference_ = pres_ref;
  }
  for (int isd = 0; isd < NumSurfDomains; isd++) {
    SurDomain1D_List[isd]->TemperatureReference_ = temp_ref;
    SurDomain1D_List[isd]->PressureReference_ = pres_ref;
  }
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
SimpleDiffusionLayout::SimpleDiffusionLayout(ProblemStatement *psInput_ptr) :
  DomainLayout(psInput_ptr), ProbNum_(0)
{
}
//=====================================================================================================================
SimpleDiffusionLayout::SimpleDiffusionLayout(int probNum,  ProblemStatement *psInput_ptr) :
  DomainLayout(psInput_ptr), 
  ProbNum_(probNum)
{
  if (probNum == 1) {
    /*
     *  Layout a simple diffusion example with trivial inlet and exit conditions
     *  This is done in a child object
     */
    InitializeDomainPicture();
  } else {
    cerr << "DomainLayout constructor: Unknown problem # " << probNum << endl;
    exit(-1);
  }

}
//=====================================================================================================================
SimpleDiffusionLayout::~SimpleDiffusionLayout() {
}
//=====================================================================================================================
SimpleDiffusionLayout::SimpleDiffusionLayout(const SimpleDiffusionLayout &r) :
  DomainLayout(), ProbNum_(0)
{
  *this = r;
}
//=====================================================================================================================
SimpleDiffusionLayout &
SimpleDiffusionLayout::operator=(const SimpleDiffusionLayout &r)
{
  if (this == &r) {
    return *this;
  }
  DomainLayout::operator=(r);

  ProbNum_ = r.ProbNum_;

  return *this;
}
//=====================================================================================================================
// Allocate the domain structure
void
SimpleDiffusionLayout::malloc_domains()
{
  int nn = psInput_ptr_->initDefaultNumCVsPerDomain_; 
  int numNodes = nn;
  double startZ = 0.0;
  double endZ = 1.0;
  BulkDomainDescription *bdd = new BDT_SimpleDiff(this, 0);

  addBulkDomainToRightEnd(bdd, numNodes, startZ, endZ);

  SurfDomainDescription *sddL = new SDT_Dirichlet(this, 1.0);
  addSurfDomainToLeftEnd(sddL, bdd);

  SurfDomainDescription *sddR = new SDT_Dirichlet(this, 0.0);
  addSurfDomainToRightEnd(sddR, bdd);

}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
SimpleTimeDependentDiffusionLayout::SimpleTimeDependentDiffusionLayout(ProblemStatement *psInput_ptr) :
  DomainLayout(psInput_ptr),
  ProbNum_(0)
{
  /*
   *  Layout a simple diffusion example with trivial inlet and exit conditions
   *  This is done in a child object
   */
  InitializeDomainPicture();
}
//======================================================================================================================
SimpleTimeDependentDiffusionLayout::SimpleTimeDependentDiffusionLayout(int probNum, ProblemStatement *psInput_ptr) :
  DomainLayout(psInput_ptr), ProbNum_(probNum)
{
  if (probNum == 1) {
    /*
     *  Layout a simple diffusion example with trivial inlet and exit conditions
     *  This is done in a child object
     */
    InitializeDomainPicture();
  } else {
    cerr << "DomainLayout constructor: Unknown problem # " << probNum << endl;
    exit(-1);
  }

}
//===========================================================================
SimpleTimeDependentDiffusionLayout::~SimpleTimeDependentDiffusionLayout() {
}
//===========================================================================
SimpleTimeDependentDiffusionLayout::SimpleTimeDependentDiffusionLayout(const SimpleTimeDependentDiffusionLayout &r) :
  DomainLayout(), ProbNum_(0)
{
  *this = r;
}
//===========================================================================
SimpleTimeDependentDiffusionLayout &
SimpleTimeDependentDiffusionLayout::operator=(const SimpleTimeDependentDiffusionLayout &r)
{
  if (this == &r) {
    return *this;
  }
  DomainLayout::operator=(r);

  ProbNum_ = r.ProbNum_;

  return *this;
}
//===========================================================================
// Allocate the domain structure
void
SimpleTimeDependentDiffusionLayout::malloc_domains()
{
  int numNodes = 10;
  double startZ = 0.0;
  double endZ = 1.0;
  BulkDomainDescription *bdd = new BDT_SimpleTDDiff(this, 0);

  addBulkDomainToRightEnd(bdd, numNodes, startZ, endZ);

  SurfDomainDescription *sddL = new SDT_Dirichlet(this, 1.0);
  addSurfDomainToLeftEnd(sddL, bdd);

  SurfDomainDescription *sddR = new SDT_Dirichlet(this, 0.0);
  addSurfDomainToRightEnd(sddR, bdd);

}
//===========================================================================
}
