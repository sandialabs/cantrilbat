/**
 * @file BaseEntry.cpp
 *   Definitions for the base object, BaseEntry
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BaseEntry.h"
#include "mdp_allo.h"
#include "BI_Dependency.h"
#include "mdp_stringUtils.h"


#ifdef WIN32
#define SNPRINTF _snprintf
#else
//! Macro to get the correct name 
#define SNPRINTF snprintf
#endif

using namespace mdpUtil;

namespace BEInput
{

/*
 *
 * BaseEntry static variables
 */
int BaseEntry::s_SkipUnknownEntries = 1;

/*
 * LineEntry static variables
 */
bool BaseEntry::s_PrintProcessedLine = true;

/*
 *  BaseEntry Constructor
 *
 * Input
 *   entryName = Name of the Base entry (required parameter)
 */
BaseEntry::BaseEntry(const char* entryName, int numTR, int timesRequiredType) :
    m_numTimesRequired(numTR),
    m_TimesRequiredType(timesRequiredType),
    m_numTimesProcessed(0),
    EntryName(entryName),
    m_multiContribIndex(0),
    NumEntryDependencies(0),
    EntryDepList(0)
{
}

/*
 *
 * BaseEntry(const BaseEntry &b):
 *
 * Copy constructor.
 * Notes -> we make sure to do a deep copy of the Dependency list so
 *          that frees will not be a problem.
 */
BaseEntry::BaseEntry(const BaseEntry& b) :
    m_numTimesRequired(b.m_numTimesRequired),
    m_TimesRequiredType(b.m_TimesRequiredType),
    m_numTimesProcessed(0),
    EntryName(b.EntryName),
    m_multiContribIndex(b.m_multiContribIndex),
    NumEntryDependencies(b.NumEntryDependencies),
    EntryDepList(0)
{
    EntryDepList = (BI_Dependency**) mdp_alloc_ptr_1(NumEntryDependencies);
    for (int i = 0; i < NumEntryDependencies; i++) {
        EntryDepList[i] = b.EntryDepList[i]->duplicateMyself();
    }
}

/*
 *  BaseEntry& operator=(const BaseEntry& b) :
 *
 * Assignment operator
 *  Notes -> we make sure to do a deep copy so that we won't have
 *           problems in the destructor.
 */
BaseEntry& BaseEntry::operator=(const BaseEntry& b)
{
    if (this != &b) {
        m_numTimesRequired  = b.m_numTimesRequired;
        m_TimesRequiredType = b.m_TimesRequiredType;
        m_numTimesProcessed = 0;
        EntryName         = b.EntryName;
        m_multiContribIndex = b.m_multiContribIndex;
        if (EntryDepList) {
            for (int i = 0; i < NumEntryDependencies; i++) {
                delete EntryDepList[i];
            }
            mdp_safe_free((void**) &EntryDepList);
        }
        NumEntryDependencies = b.NumEntryDependencies;
        EntryDepList = (BI_Dependency**) mdp_alloc_ptr_1(NumEntryDependencies);
        for (int i = 0; i < NumEntryDependencies; i++) {
            EntryDepList[i] = b.EntryDepList[i]->duplicateMyself();
        }
    }
    return *this;
}

/*
 *
 * duplMyselfAsBaseEntry() (virtual) :
 *
 * Create a duplication function that uses the copy constructor.
 * This function, in contrast to the copy constructor and the
 * operator= function can be virtual. Therefore we can copy lists
 * of BaseEntry's while retaining derived type information by
 * overriding this function in the derived classes.
 */
BaseEntry* BaseEntry::duplMyselfAsBaseEntry() const
{
    BaseEntry* newBE = new BaseEntry(*this);
    return newBE;
}

/*
 *
 * ~BaseEntry() (virtual function) :
 */
BaseEntry::~BaseEntry(void)
{
    /*
     * Dependency objects are owned by the Entries. Therefore,
     * when we delete the entries, we must delete the
     * dependency objects and the pointer list used to hold them.
     */
    if (EntryDepList) {
        for (int i = 0; i < NumEntryDependencies; i++) {
            delete EntryDepList[i];
        }
    }
    mdp_safe_free((void**) &EntryDepList);
}

/*
 * zeroTimesCounter():
 *
 *  This call zeroes the line count information out.
 */
void BaseEntry::zeroTimesCounter()
{
    m_numTimesProcessed = 0;
}

void BaseEntry::set_NumTimesRequired(int ntr)
{
    m_numTimesRequired = ntr;
}

void BaseEntry::set_TimesRequiredType(int timesRequiredType)
{
    m_TimesRequiredType = timesRequiredType;
}

int BaseEntry::get_NumTimesProcessed() const
{
    return m_numTimesProcessed;
}

std::string BaseEntry::keyname() const
{
    return std::string(EntryName.orig_str);
}

// Return the multiContribIndex value
int BaseEntry::multiContribIndex() const
{
    return m_multiContribIndex;
}

void BaseEntry::set_multiContribIndex(int contribValue)
{
    m_multiContribIndex = contribValue;
}
/*
 * Declare a dependency for this entry on another entry.
 * And declare a type.
 */
void BaseEntry::declareDependency(BI_Dependency* bi_dep)
{
    if (bi_dep) {
        /*
         *  Get the type of dependency
         */
        BIDSR_TYPE srtype = bi_dep->ServiceRequestType();
        const BaseEntry* be = bi_dep->TargetBaseEntry();
        if (!be) {
            throw BI_InputError("BaseEntry::declareDependency",
                                "Required BaseEntry is null");
        }
        /*
         * Ask the object if it can service the dependency
         * If it can, we are good to go.
         * We add this dependency into the internal list
         * of dependencies for this object.
         */
        if (be->DepCanService(srtype)) {
            mdp_realloc_ptr_1((void***) &EntryDepList,
                              NumEntryDependencies + 1,
                              NumEntryDependencies);
            EntryDepList[NumEntryDependencies] = bi_dep;
            NumEntryDependencies++;
        } else {
            throw BI_InputError("BaseEntry::declareDependency",
                                "object can't service request type");
        }
    } else {
        BI_InputError("BaseEntry::declareDependency",
                      "bi_dep is NULL");
    }
}

/*
 *
 * declareSimpleOrderDependency():
 *
 * Declare a simple order dependency.
 */
void BaseEntry::declareSimpleOrderDependency(BaseEntry* be_req)
{
    if (be_req) {
        /*
         *
         */
        BIDT_TYPE dtype =  BIDT_ENTRYPROCESSED;
        BI_Dependency* bi_dep = new BI_Dependency(be_req, dtype,
                                                  BIDRT_PT_ERROR);
        declareDependency(bi_dep);
    } else {
        BI_InputError("BaseEntry::declareSimpleOrderDependency",
                      "be_req is NULL");
    }
}

/*
 *
 * DepCanService():
 *
 *  This call checks to see whether
 *
 * @param BIDSR_value  Dependency service request type
 *      -> right now this is figured out from the previous two ints.
 *         It is not part of the interface.
 */
bool BaseEntry::DepCanService(BIDSR_TYPE BIDSR_value) const
{
    if (BIDSR_value == BIDSR_SIMPLEORDERING) {
        return true;
    }
    return false;
}

/*
 * ansDepCheck():
 *
 *  This routine answers a request from another BaseEntry object as to
 *  whether this object has fulfilled its requirements in the
 *  input deck at the current point in processing the input deck.
 *
 *  If the entry is not required, then it will return true if it has
 *  been called at all. If it has not been called, it will return false.
 *  If the entry is required, it will return true if the number
 *  of calls equals or exceeds the required number of calls.
 */
bool BaseEntry::ansDepCheck() const
{
    if (!m_numTimesRequired) {
        if (m_numTimesProcessed) {
            return true;
        }
    } else {
        if (m_numTimesProcessed >= m_numTimesRequired) {
            return true;
        }
    }
    return false;
}

/*
 *
 * checkDepOneInt(): (virtual)
 *
 * This routine answers the same dependency check as ansDepCheck
 * and it returns a single int value depending upon the
 * value input in the input file.
 */
bool BaseEntry::ansDepCheckOneInt(int& returnInt) const
{
    throw BI_InputError("BaseEntry::ansDepCheckOneInt",
                        "object can't service dependency type");
}

// Return an address from which an application can read the
// data, if it knows how to convert the memory address.
/*
 * This address is a const void * since we don't know what
 * the form of the data is.
 */
const void* BaseEntry::currentValueAsVoidP() const
{
    return static_cast<const void*>(0);
}

bool BaseEntry::get_printProcessedLine()
{
    return s_PrintProcessedLine;
}
void BaseEntry::set_printProcessedLine(bool bv)
{
    s_PrintProcessedLine = bv;
}

int BaseEntry::get_SkipUnknownEntries()
{
    return s_SkipUnknownEntries;
}
void BaseEntry::set_SkipUnknownEntries(int bv)
{
    s_SkipUnknownEntries = bv;
}
}
