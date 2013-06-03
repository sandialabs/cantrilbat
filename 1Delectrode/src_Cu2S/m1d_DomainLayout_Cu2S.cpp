/**
 * @file m1d_DomainLayout_Cu2S.cpp
 *
 */

/*
 *  $Id: m1d_DomainLayout_Cu2S.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */
#include "m1d_defs.h"
#include "m1d_DomainLayout.h"
#include "m1d_BulkDomainTypes.h"
#include "m1d_SurfDomainTypes.h"

#include "m1d_DomainLayout_Cu2S.h"
#include "m1d_BDT_Cu2S.h"
#include "m1d_SDT_Mixed_Cu2S.h"

#include "m1d_exception.h"

#include <iostream>

using namespace std;

namespace m1d
{

//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
DomainLayout_Cu2S::DomainLayout_Cu2S() :
  DomainLayout(), ProbNum_(0)
{
  /*
   *  Layout a simple diffusion example with trivial inlet and exit conditions
   *  This is done in a child object
   */
  InitializeDomainPicture();
}
//===========================================================================
DomainLayout_Cu2S::DomainLayout_Cu2S(int probNum, ProblemStatement *psInput_ptr) :
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
DomainLayout_Cu2S::~DomainLayout_Cu2S()
{
}
//===========================================================================
DomainLayout_Cu2S::DomainLayout_Cu2S(const DomainLayout_Cu2S &r) :
  DomainLayout(), ProbNum_(0)
{
  *this = r;
}
//===========================================================================
DomainLayout_Cu2S &
DomainLayout_Cu2S::operator=(const DomainLayout_Cu2S &r)
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
DomainLayout_Cu2S::malloc_domains()
{
  int numNodes = 10;
  double startZ = 0.0;
  double endZ = 1.0E-6;
  BulkDomainDescription *bdd = new BDT_Cu2S(this, 0);

  addBulkDomainToRightEnd(bdd, numNodes, startZ, endZ);

  VarType sp0(Concentration_Species, 0);
  EqnType sp0_cons(Species_Conservation, 0);

  // SDT_Dirichlet*dirLeft = new SDT_Dirichlet(this);
  // dirLeft->addDirichletCondition(sp0_cons, sp0, 1.0);
  SDT_Mixed_Cu2S * dirLeft = new SDT_Mixed_Cu2S(this, 1);
  SurfDomainDescription *sddL = dirLeft;
  addSurfDomainToLeftEnd(sddL, bdd);

  SDT_Mixed_Cu2S * dirRight = new SDT_Mixed_Cu2S(this, 0);
  //dirRight->addDirichletCondition(sp0, 0.0);
  SurfDomainDescription *sddR = dirRight;
  addSurfDomainToRightEnd(sddR, bdd);

}
//===========================================================================
}
