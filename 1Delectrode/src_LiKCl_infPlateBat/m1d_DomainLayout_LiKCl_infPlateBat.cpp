/**
 * @file m1d_DomainLayout_LiKCl_infPlateBat.cpp
 *
 */

#include "m1d_defs.h"
#include "m1d_DomainLayout.h"
#include "m1d_BulkDomainTypes.h"
#include "m1d_SDD_Mixed.h"

#include "m1d_DomainLayout_LiKCl_infPlateBat.h"
#include "m1d_BDD_porousLiKCl.h"
#include "m1d_SDD_FlatAnode.h"
#include "m1d_SDD_FlatCathode.h"

#include "m1d_exception.h"

#include <iostream>

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace std;

namespace m1d
{

  //=====================================================================================================================
  //=====================================================================================================================
  //=====================================================================================================================
  DomainLayout_LiKCl_infPlateBat::DomainLayout_LiKCl_infPlateBat(ProblemStatement* psInput_ptr) :
    DomainLayout(psInput_ptr), ProbNum_(0), pscInput_ptr_(0)
  {
     pscInput_ptr_ = dynamic_cast<ProblemStatementCell *>(psInput_ptr);
     if (!pscInput_ptr_) {
      Zuzax::ZuzaxError("DomainLayout_LiKCl_infPorousBat::DomainLayout_LiKCl_infPorousBat()",
                          "Bad dynamic cast");
     } 
     InitializeDomainPicture();
  }
  //===========================================================================
  DomainLayout_LiKCl_infPlateBat::DomainLayout_LiKCl_infPlateBat(int probNum, ProblemStatement *psInput_ptr) :
    DomainLayout(psInput_ptr), ProbNum_(probNum), pscInput_ptr_(0)
  {
    pscInput_ptr_ = dynamic_cast<ProblemStatementCell *>(psInput_ptr);
    if (!pscInput_ptr_) {
      Zuzax::ZuzaxError("DomainLayout_LiKCl_infPorousBat::DomainLayout_LiKCl_infPorousBat()",
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
  DomainLayout_LiKCl_infPlateBat::~DomainLayout_LiKCl_infPlateBat()
  {
  }
  //====================================================================================================================
  DomainLayout_LiKCl_infPlateBat::DomainLayout_LiKCl_infPlateBat(const DomainLayout_LiKCl_infPlateBat &r) :
    DomainLayout(), ProbNum_(0)
  {
    *this = r;
  }
  //====================================================================================================================
  DomainLayout_LiKCl_infPlateBat &
  DomainLayout_LiKCl_infPlateBat::operator=(const DomainLayout_LiKCl_infPlateBat &r)
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
  DomainLayout_LiKCl_infPlateBat::malloc_domains()
  {

    /*
     *  Here we lay out the size of the domain
     */
    double startZ = 0.0;
    //    double endZ = 0.4826E-3; //actual separator thickness
    double endZ = PSinput.separatorThickness_; 
  
    BulkDomainDescription *bdd = new BDD_porousLiKCl(this);

    int numNodesS = pscInput_ptr_->initDefaultNumCVsSeparator_;
    addBulkDomainToRightEnd(bdd, numNodesS, startZ, endZ);

    ///VarType sp0(Concentration_species, 0);
    // EqnType sp0_cons(Species_conservation, 0);


    SDD_FlatAnode * dirLeft = new SDD_FlatAnode(this, 1);
    SurfDomainDescription *sddL = dirLeft;
    addSurfDomainToLeftEnd(sddL, bdd);

    //Here if we want to specify a current profile we can do so
    //with an alternate ProbNum_ and an alternate SDT_FlatCathode()

    SDD_FlatCathode * dirRight = new SDD_FlatCathode(this, 1);

    SurfDomainDescription *sddR = dirRight;
    addSurfDomainToRightEnd(sddR, bdd);

  }
  //===========================================================================
}
