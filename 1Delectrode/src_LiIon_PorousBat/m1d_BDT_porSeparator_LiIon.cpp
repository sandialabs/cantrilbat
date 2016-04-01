/**
 * @file m1d_BDT_porousSeparator_LiIon.cpp
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_BDT_porSeparator_LiIon.h"
#include "m1d_porousLiIon_Separator_dom1D.h"
#include "m1d_ProblemStatementCell.h"
#include "LiIon_PorousBat.h"

using namespace std;

namespace m1d
{

//==================================================================================================================================
BDT_porSeparator_LiIon::BDT_porSeparator_LiIon(DomainLayout* dl_ptr, std::string domainName) :
     BDD_porousFlow(dl_ptr, "separator", domainName)
{
    //int eqnIndex = 0;
    IsAlgebraic_NE.resize(6,0);
    IsArithmeticScaled_NE.resize(6,0);
}
//===================================================================================================================================
BDT_porSeparator_LiIon::BDT_porSeparator_LiIon(const BDT_porSeparator_LiIon& r) :
    BDD_porousFlow(r.DL_ptr_)
{
    *this = r;
}
//======================================================================================================================================
BDT_porSeparator_LiIon::~BDT_porSeparator_LiIon()
{
}
//======================================================================================================================================
BDT_porSeparator_LiIon&
BDT_porSeparator_LiIon::operator=(const BDT_porSeparator_LiIon& r)
{
    if (this == &r) {
        return *this;
    }

    BDD_porousFlow::operator=(r);

    return *this;
}
//======================================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
BulkDomain1D*
BDT_porSeparator_LiIon::mallocDomain1D()
{
    BulkDomainPtr_ = new porousLiIon_Separator_dom1D(this);
    return BulkDomainPtr_;
}
//======================================================================================================================================
} /* End of Namespace */
//==================================================================

