/**
 * @file m1d_DomainLayout_Cu2S.cpp
 *
 */

#include "m1d_defs.h"
#include "m1d_DomainLayout_Cu2S.h"
#include "m1d_BDD_Cu2S.h"
#include "m1d_SDD_Mixed_Cu2S.h"
#include "m1d_DomainLayout.h"

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
  BulkDomainDescription *bdd = new BDD_Cu2S(this, 0);

  addBulkDomainToRightEnd(bdd, numNodes, startZ, endZ);

  VarType sp0(Concentration_Species, 0);
  EqnType sp0_cons(Species_Conservation, 0);

  // SDD_Dirichlet*dirLeft = new SDD_Dirichlet(this);
  // dirLeft->addDirichletCondition(sp0_cons, sp0, 1.0);
  SDD_Mixed_Cu2S * dirLeft = new SDD_Mixed_Cu2S(this, 1);
  SurfDomainDescription *sddL = dirLeft;
  addSurfDomainToLeftEnd(sddL, bdd);

  SDD_Mixed_Cu2S * dirRight = new SDD_Mixed_Cu2S(this, 0);
  //dirRight->addDirichletCondition(sp0, 0.0);
  SurfDomainDescription *sddR = dirRight;
  addSurfDomainToRightEnd(sddR, bdd);

}
//===========================================================================
}
