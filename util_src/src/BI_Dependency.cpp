/**
 * @file BI_Dependency.cpp
 *    Definitions for the BI_Dependency class
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BaseEntry.h"
#include "BI_Dependency.h"

#include <new>
namespace BEInput {
  /*
   * BI_Dependency():
   *
   * Default constructor
   */
  BI_Dependency::BI_Dependency(const BaseEntry *be, 
			       BIDT_TYPE BIDT_type,
			       BIDRT_TYPE BIDRT_type) :		     
    RequiredEntry(be),
    m_DT_type(BIDT_type),
    m_RT_type(BIDRT_type),
    m_SR_type(BIDSR_SIMPLEORDERING),
    CardMax( 1000000),
    CardMin(-1000000)
  {
    if (m_DT_type == BIDT_ONEINT) {
      m_SR_type = BIDSR_ONEINT;
    }
    if (m_DT_type == BIDT_INTMAXMIN) {
      m_SR_type = BIDSR_ONEINT;
    }
    if (!be) {
      throw BI_InputError("BI_Dependency", "target BaseEntry is Null");
    }
  }

  /*
   *
   *   BI_Dependency(const BI_Dependency&):
   *
   *  copy constructor
   */
  BI_Dependency::BI_Dependency(const BI_Dependency& b) :
    RequiredEntry(b.RequiredEntry),
    m_DT_type(b.m_DT_type),
    m_RT_type(b.m_RT_type),
    m_SR_type(b.m_SR_type),
    CardMax(b.CardMax),
    CardMin(b.CardMin)
  {
  }

  /*
   *   BI_Dependency& operator=(const BI_Dependency&):
   *
   * assigntment operator
   */
  BI_Dependency& BI_Dependency::operator=(const BI_Dependency&b) {
    if (&b != this) {
      RequiredEntry = b.RequiredEntry;
      m_DT_type     = b.m_DT_type;
      m_RT_type     = b.m_RT_type;
      m_SR_type     = b.m_SR_type;
      CardMax       = b.CardMax;
      CardMin       = b.CardMin;
    }
    return *this;
  }
  /*
   * ~BI_Dependency():
   *
   *  Destructor
   */
  BI_Dependency::~BI_Dependency() {
  }

  /*
   * BI_Dependency* duplicateMyself();  (virtual)
   *
   * Duplication function
   */
  BI_Dependency* BI_Dependency::duplicateMyself() const {
    BI_Dependency* newDep = new BI_Dependency(*this);
    return newDep;
  }

  /*
   * BI_Dependency::checkDependencySatisfied():
   *
   */
  bool BI_Dependency::checkDependencySatisfied() const {
    bool retn = RequiredEntry->ansDepCheck();
    return retn;
  }

  /*
   * BI_Dependency::checkDepOneInt():
   *
   */
  bool BI_Dependency::checkDepOneInt(int &retnInt) const {
    throw BI_InputError("BI_Dependency::checkDepOneInt",
			"not handled");
  }

  /*
   * BI_Dependency::checkCardIntMaxMin():
   *
   */
  bool BI_Dependency::checkCardIntMaxMin(int cardInt) {
    if (cardInt > CardMax) return false;
    if (cardInt < CardMin) return false;
    return true;
  }

  /*
   * BI_Dependency::setCardIntMaxMin():
   *
   */
  void BI_Dependency::setCardIntMaxMin(int cMax, int cMin) {
    CardMax = cMax;
    CardMin = cMin;
  }

  BIDRT_TYPE BI_Dependency::ResultType() const {
    return m_RT_type;
  }

  BIDSR_TYPE BI_Dependency::ServiceRequestType() const {
    return m_SR_type;
  }

  const BaseEntry * BI_Dependency::TargetBaseEntry() const {
    return RequiredEntry;
  }
}
