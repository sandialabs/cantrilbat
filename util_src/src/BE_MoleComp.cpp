/**
 * @file BE_MoleComp.cpp
 *  Definitions for the BlockEntry of a set of doubles that fills up a vector
 *   of mole fractions
 *  (see \ref blockentryModule and class
 *  \link BEInput::BE_MoleComp BE_MoleComp\endlink).
 */
/*
 * $Author: hkmoffa $
 * $Revision: 5 $
 * $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_MoleComp.h"

using namespace BEInput;
namespace BEInput
{
/*
 * BE_MoleComp Constructor:
 *
 *   This sets up the line entry special case.
 *   We make sure to call the base class constructor here to do
 *   much of the initialization.
 */
BE_MoleComp::BE_MoleComp(const char* blockName, double** hndlAAA,
                         int numTimesRequired,
                         char** charList, int listLength,
                         int constructLE,
                         const char* varName,
                         BlockEntry* parentBlock_input) :
    BE_StrDbl(blockName, hndlAAA, numTimesRequired, 0,
              charList,listLength, constructLE,
              varName, parentBlock_input)
{
    /*
     * We set the default to 0. This means that the mole fraction
     * vector is zeroed just before being filled.
     */
    set_default(0.0);
    set_limits(1.0, 0.0);
}

/*
 * BE_MoleComp(const BE_MoleComp&)
 *
 * Copy constructor
 */
BE_MoleComp::BE_MoleComp(const BE_MoleComp& b) :
    BE_StrDbl(b)
{
}

/*
 *  BE_MoleComp& operator=(const BE_MoleComp&);
 *
 *  Assignment operator
 */
BE_MoleComp& BE_MoleComp::operator=(const BE_MoleComp& b)
{
    if (&b != this) {
        BE_StrDbl::operator=(b);
    }
    return *this;
}

/*
 * BlockEntry* duplMyselfAsBlockEntry();
 *
 *  Duplicate myself in a list of base class objects
 */
BlockEntry* BE_MoleComp::duplMyselfAsBlockEntry() const
{
    BE_MoleComp* newBE = new BE_MoleComp(*this);
    return (BlockEntry*) newBE;
}

/*
 * BE_MoleComp destructor: (virtual function)
 *
 * We malloced memory here, so we must explicitly call free.
 */
BE_MoleComp::~BE_MoleComp()
{
#ifdef DEBUG_DESTRUCTOR
    printf("~BE_MoleComp called for %s\n", BlockName.orig_str);
#endif
}

/*
 *  Wrapup() (virtual function)
 *
 *  This function is called when the end block line for the
 *  current block is read. Cleanup is done, and debugging printouts
 *  as well.
 *  We require the mole fractions to sum to almost 1.0. Small errors
 *  are gotten rid of by normalization. Large errors produce
 *  errors.
 */
void BE_MoleComp::Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    /*
     * Normalize
     */
    double sum = 0.0;
    double* AddrVal = *HndlDblVec;
    for (int i = 0; i < ListLength; i++) {
        sum += AddrVal[i];
    }
    if (sum <= 0.98 || sum >= 1.02) {
        char st[20];
        sprintf(st, "%g", sum);
        std::string err = "Incorrect total sum of input mole fractions: sum = " ;
        err += st;
        throw BI_InputError("BE_MoleComp::Wrapup", err);
    } else {
        for (int i = 0; i < ListLength; i++) {
            AddrVal[i] /= sum;
        }
    }

    BE_StrDbl::Wrapup(ifp_input, blockArgTok);
}
}
