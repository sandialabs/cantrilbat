/**
 * @file BaseEntry.h
 *   Declarations for the base object, BaseEntry
 */
/* $Author: hkmoffa $
 * $Revision: 508 $
 * $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#ifndef BASEENTRY_H
#define BASEENTRY_H

/**
 *  LONG_PTR should hold the value of the pointer, and be able
 *     to do arithmetic with the pointer. Note, since pointers are
 *     basically unsigned ints, this means that long long should be
 *     used when possible.
 *     In other words, if pointers are 4 bytes, then a 4 byte signed
 *     int may not be enough to hold all values, in all cases. You
 *     should go to 8 bytes.
 */
#ifndef LONG_PTR
#ifdef WIN32
#define LONG_PTR long long
#else
#define LONG_PTR long long
#endif
#endif

#include <cstdlib>
#include <cstdio>
#include "tok_input_util.h"
#include "BI_InputError.h"
#include "BI_Dependency.h"

namespace BEInput {
/*
 * Forward reference needed here.
 */
class BI_Dependency;

//! Typedef for the token structure used so frequently
/*!
 * Bringing it into the BEInput namespace as it
 * is used in the interface to public functions
 */
typedef struct TKInput::TOKEN TK_TOKEN;

//! Convert a double into a c++ string
/*
 *  This routine doesn't assume a formatting. You
   *  must supply the formatting
   *
   * @param x double to be converted
   * @param fmt   Format to be used (printf style)
   */
std::string fp2str(const double x);

//!   This is the base class for all LineEntry and BlockEntry objects.
//!   At this level dependencies between objects are handled 
/*!
 *   The member EntryName contains the keyline entry name of blocks and lines
 *    
 * @ingroup blockentryModule
 */
class BaseEntry {

public:

  //! Base Constructor
  /*!
   *   @param entryName  keyline Name of the Base entry (required parameter)
   *   @param numTR      Number of times this card must be found 
   *                     in the input deck.
   *   @param timesRequiredType Sets the value of m_TimesRequiredType in the
   *                            object. This determines the type of requirement
   *                            for the number of times processed.
   */
  explicit BaseEntry(const char *entryName, int numTR = 0, int timesRequiredType = 0);

  //! Copy Constructor
  /*!
   * @param right object to be copied
   */
  BaseEntry(const BaseEntry& right);

  //! Assignment operator
  /*!
   * @param right object to be copied
   */
  BaseEntry& operator=(const BaseEntry& right);

  //! Destructor
  virtual ~BaseEntry();

  //! Duplicator function.
  /*!
   * This function duplicates the entry and returns a pointer
   * to a BaseEntry. 
   */
  virtual BaseEntry* duplMyselfAsBaseEntry() const;
    
  //! This call zeroes the NumTimesProcessed information out. 
  virtual void zeroTimesCounter();

  //! Set the number of times this card is required in the
  //! input deck
  /*!
   *  This can be used to override the entry in the constructor
   *
   * @param ntr Int value for the number of times this 
   *            card is required
   */
  void set_NumTimesRequired(int ntr);

  //! Set the times required type field which determines how
  //! the NumTimesRequired field will be handled
  /*!
   *   There are four types
   *          0 = exact match for the number of times required
   *          1 = greater than match:  The number of times processed must
   *                    be greater than or equal to m_numTimesRequired.
   *          2 = less than or equal match: The number of times processed
   *                    must be less than or equal to m_numTimesRequired.
   *                    It may be zero.
   *          3 = less than or equal match: The number of times processed
   *                    must be less than or equal to m_numTimesRequired.
   *                    It must be at least one.
   */
  void set_TimesRequiredType(int timesRequiredType);

  //! return the number of times processed int
  int get_NumTimesProcessed() const;

  //! return the name of the object
  std::string keyname() const;
    
  //! Declare a dependency for this entry on another entry.
  /*!
   * The dependency has a type and parameters. But this function only
   * employs the base dependency class.
   *
   * Ask the Queried BaseEntry object if it can service the dependency. If it can, we are good to go.
   * We add this dependency into the internal list of dependencies for this object.
   * The BaseEntry owns the dependency object and will try to delete it during the
   * destructor.
   *
   * @param bi_dep Pointer to the dependency object that describes the dependency.
   */
  void declareDependency(BI_Dependency *bi_dep);

  //! Declare the arugment BaseEntry must have occurred before this BaseEntry
  /*!
   *  This dependency can be used to force ordering within the
   *  blocks and line entries of the input deck
   *
   * @param be_req Pointer to an existing BaseEntry.
   */
  void declareSimpleOrderDependency(BaseEntry *be_req);

  //! Return the multiContribIndex value
  int multiContribIndex() const;

  void set_multiContribIndex(int contribIndex);

  //! Query to this object, asking if this object can service
  //! a particular dependency service request type.
  /*!
   * @param BIDSR_value  Dependency service request type
   *      -> right now this is figured out from the previous two ints.
   *         It is not part of the interface.
   */
  virtual bool DepCanService(BIDSR_TYPE BIDSR_value) const;

  //!  This routine answers a request from another BaseEntry object as to
  //!  whether this object has fulfilled its requirements in the
  //!  input deck at the current point in processing the input deck. 
  /*!
   *  If the entry is not required, then it will return true if it has
   *  been called at all. If it has not been called, it will return false.
   *  If the entry is required, it will return true if the number
   *  of calls equals or exceeds the required number of calls.
   */
  bool ansDepCheck() const;
  
  //! This routine answer the same dependency check as ansDepCheck(),
  //! and it returns a single int value, which is the
  //! value input in the input file pertaining to this card.
  /*!
   *  @param returnInt Return value of the int.
   *                   This int is the natural int value
   *                   Created by the derived class of this int.           
   *  @return
   *     This returns the same value as the function
   *     ansDepCheck() does.
   */
  virtual bool ansDepCheckOneInt(int &returnInt) const;

  //! Return an address from which an application can read the
  //! data, if it knows how to convert the memory address.
  /*!
   * This address is a const void * since we don't know what
   * the form of the data is.
   * For the base case, we return the NULL pointer.
   *
   * @note Can not make this pure virtual yet, because
   *       this functions hasn't been implemented in all objects yet.
   */
  virtual const void *currentValueAsVoidP() const;

  //! get the Boolean indicating whether a line should be echoed to output
  //! if it is found and processed
  static bool get_printProcessedLine();

  //! Set the Boolean indicating whether a line should be echoed to output
  //! if it is found and processed
  /*!
   * @param bv value to be set
   */
  static void set_printProcessedLine(bool bv);

  //! get the Boolean indicating  whether an error is thrown if an unknown
  //!  block or line entry is encountered in the input.
  static bool get_SkipUnknownEntries();

  //! Set the Boolean indicating  whether an error is thrown if an unknown
  //! block or line entry is encountered in the input.
  /*!
   * @param bv value to be set
   */
  static void set_SkipUnknownEntries(bool bv);

protected:
  //!  Number of times this object is required to be found in the
  //!  input deck for a proper specification of the problem definition
  /*!
   * The default value is zero.
   */
  int m_numTimesRequired;

  //!  Type of the of the numTimesRequired condition
  /*!
   *   There are three times
   *          0 = exact match for the number of times required
   *          1 = greater than match:  The number of times processed must
   *                    be greater than or equal to m_numTimesRequired.
   *          2 = less than or equal match: The number of times processed
   *                    must be less than or equal to m_numTimesRequired.
   *                    It may be zero.
   *          3 = less than or equal match: The number of times processed
   *                    must be less than or equal to m_numTimesRequired.
   *                    It must be at least one.
   *   The default is type 0;
   */
  int m_TimesRequiredType;
  
  //! Number of times this object has been called for the current parent
  //! block.
  int m_numTimesProcessed;
  
  //! This is the keyline name of the BaseEntry
  TK_TOKEN EntryName;

  //! Index of the Base entry, for the case where there
  //! can be multiple BlockEntries
  /*!
   *  This is almost always 0. When there are multiblocks
   *  each multiblock will have a unique id, but
   *  the same EntryName.
   */
  int m_multiContribIndex;

  //! Number of Dependencies for this entry.
  /*!
   * When this entry is utilized, there may be dependencies on 
   *    other entries having already been processed. This field
   *    yields the number of these dependencies.
   *    It is an error for an entry to be called without first
   *    having all of its dependencies fulfilled.
   */
  int NumEntryDependencies;

  //!  List of dependencies. Note, this object is the owner of this list
  //!  of dependencies. It will free them when this object goes away.
  /*!
   * The length of the list is equal to NumEntryDependencies
   */
  BI_Dependency** EntryDepList;

  //! Boolean indicating whether a line should be echoed to output
  //! if it is found and processed
  /*!
   * The default value for this is true
   */
  static bool s_PrintProcessedLine;

  //! Boolean indicating  whether an error is thrown if an unknown
  //! block or line entry is encountered in the input.
  /*!
   * The default is false. -> throw an error if unknown entries
   * are encountered.
   */
  static bool s_SkipUnknownEntries;

};
}
/***************************************************************************/
#endif
