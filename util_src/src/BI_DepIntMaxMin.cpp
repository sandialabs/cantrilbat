/**
 * @file BI_DepIntMaxMin.cpp
 *
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

#include "BI_DepIntMaxMin.h"
#include "mdp_allo.h"
#include "BI_Dependency.h"
#include <new>

namespace BEInput
{

/*
 * BI_DepIntMaxMin() constructor
 *
 */
BI_DepIntMaxMin::BI_DepIntMaxMin(BaseEntry* be, BIDT_TYPE BIDT_type,
                                 int maxTargetI, int minTargetI,
                                 BIDRT_TYPE BIDR_type) :
    BI_Dependency(be, BIDT_type, BIDR_type),
    m_maxValue(maxTargetI),
    m_minValue(minTargetI),
    m_RetnOneInt(MDP_INT_NOINIT)
{
    m_SR_type = BIDSR_ONEINT;
    if (BIDT_type != BIDT_ONEINT && BIDT_type != BIDT_INTMAXMIN) {
        throw BI_InputError("BI_DepIntMaxMin",
                            "BIDT_TYPE value is not appropriate");
    }
    if (m_maxValue < m_minValue) {
        throw BI_InputError("BI_DepIntMaxMin",
                            "max value is less than min value");
    }
}

/*
 *   BI_DepIntMaxMin(const BI_DepIntMaxMin&):
 *
 * copy constructor
 */
BI_DepIntMaxMin::BI_DepIntMaxMin(const BI_DepIntMaxMin& b) :
    BI_Dependency(b),
    m_maxValue(b.m_maxValue),
    m_minValue(b.m_minValue),
    m_RetnOneInt(b.m_RetnOneInt)
{
}

/*
 *   BI_DepIntMaxMin& operator=(const BI_DepIntMaxMin&):
 *
 *   assigntment operator
 */
BI_DepIntMaxMin& BI_DepIntMaxMin::operator=(const BI_DepIntMaxMin& b)
{
    if (&b != this) {
        BI_Dependency::operator=(b);
        m_maxValue = b.m_maxValue;
        m_minValue = b.m_minValue;
        m_RetnOneInt = b.m_RetnOneInt;
    }
    return *this;
}

/*
 * ~BI_DepIntMaxMin()
 *
 * destructor
 */
BI_DepIntMaxMin::~BI_DepIntMaxMin()
{
}

/*
 *  BI_Dependency* duplicateMyself
 *
 *  Duplication of object as the base class
 */
BI_Dependency* BI_DepIntMaxMin::duplicateMyself() const
{
    BI_DepIntMaxMin* newDep = new BI_DepIntMaxMin(*this);
    return (BI_Dependency*) newDep;
}

/*
 * BI_DepIntMaxMin::checkDependencySatisfied():
 *
 */
bool BI_DepIntMaxMin::checkDependencySatisfied() const
{
    bool retn = RequiredEntry->ansDepCheckOneInt(m_RetnOneInt);
    if (!retn) {
        return false;
    }
    if (m_RetnOneInt > m_maxValue) {
        return false;
    }
    if (m_RetnOneInt < m_minValue) {
        return false;
    }
    return true;
}

/*
 * BI_Dependency::checkDepOneInt():   (virtual)
 */
bool BI_DepIntMaxMin::checkDepOneInt(int& retnInt) const
{
    bool retn = RequiredEntry->ansDepCheckOneInt(m_RetnOneInt);
    retnInt = m_RetnOneInt;
    return retn;
}


}
