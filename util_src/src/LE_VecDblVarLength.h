/**
 * @file LE_VecDblVarLength.h
 *  Header for the LineEntry of a vector of  doubles
 *  (see \ref blockentryModule and class 
 *  \link BEInput::LE_VecDbl LE_VecDbl\endlink).
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

#ifndef LEVECDBLVARLENGTH_H
#define LEVECDBLVARLENGTH_H

#include "LineEntry.h"

namespace BEInput {


  //!  Keyline Entry for a variable vector of doubles.
  /*!
   *   This sets up the <b>LineEntry</b> special case for the entry of a
   *   vector of doubles. It declares interfaces for specifying the
   *   entry and for processing the entry once a match is
   *   made in the input deck. The general form of the input is
   *
   *   \f[
   *     KeyName = \mathrm{double\_1} \mathrm{double\_2} \mathrm{...}
   *   \f]
   *
   *  The <I>KeyName</I> variable  may be consist of multiple
   *  tokens. It is tokenized and compared in a case insensitive format
   *  at the start of the processing and before any matching
   *  is carried out. The rhs consists of multiple doubles.
   *  The number of doubles is a variable value, not known at 
   *  construction time.
   *  The length of the vector may be determined via post-processing
   *  by querying the number of times the card has been processed.
   *  Multiple doubles on a single card is counted as multiple times
   *  processed.
   *
   *  Additionally, the following tokens are understood and accepted in
   *  place of \f$ \mathrm{double} \f$ : <tt>DBL_MAX</tt>,  <tt>DBL_MIN</tt>,
   *  and <tt>Default</tt>. If <tt>DBL_MAX</tt> or  <tt>DBL_MIN</tt> is
   *  encountered, the value is set to the respective value as defined
   *  in the file math.h. If the value <tt>Default</tt> is encountered,
   *  the #CurrValue is set to the default value, #DefaultVal,
   *  specified for the object.
   *
   *  If the default value specified for the object 
   *  is set to <tt>NO_DEFAULT_DBL</tt>, which is the initial value for
   *  #DefaultVal, then an error is thrown when <tt>default</tt> 
   *  is encountered on the rhs.
   *
   *  When a match is found, the rhs of the card is processed into a double
   *  using the function #process_LineEntry().
   *  The double is then assigned to the #CurrValue of this object.
   *  And, it is optionally written out to an external address of a vector
   *  of doubles, #HndlDblVec,
   *  that was supplied during the construction of the object.
   *  The object may then later be queried, using either the
   *  #currentTypedValue() or #currentValueAsVoidP()
   *  functions for the value of #CurrValue
   *  during subsequent processing.
   *
   *  A maximum and minimum allowed double value may be set as a requirement.
   *
   * <H2> Example of Usage </H2>
   *
   *
   *  Example:
   *
   *   The example below sets up a required keyline entry, putting the value in
   *   the global address globInput.molarVolumes[]
   *
   * @code
   *      struct globInput {
   *         int numSpecies;
   *         double *molarVolumes;
   *      } gI;
   *
   *      BlockEntry *besmd = new BlockEntry("Fluid Properties");
   *      int reqd = 0;
   *      LE_VecDbl *v1 = new LE_VecDblVarLength("Molar Volumes",
   *                                  &(gI.molarVolumes),
   *                                    gI.numSpecies, reqd,
   *                                    "molarVolumes");
   *      v1->set_default(1.3E-5);
   *      v1->set_limits(100000., 1.0E-4);
   *      besmd->addLineEntry(v1);
   * @endcode
   *
   *  An example of the input deck entry for this one LineEntry follows.
   *  Note, the line may be enclosed in any number of nested blocks,
   *  as long as the last nested block is named "Fluid Properties".
   *  All white space differences and capitalization differences
   *  are ignored.
   *
   *  @code
   *       Start Block Fluid Properties
   *         Molar Volumes = 1.3E-5 1.3E-5 2.0E-5
   *       End Block   Fluid Properties
   *  @endcode
   *
   *  
   *  Examples of executing the code and extracting values follow.
   *
   *  The input file containing the entry may be processed via a
   *  command similar to the following
   *
   *  @code
   *     FILE *ifp = fopen("inputFile.txt", "r");
   *     besmd->read_block(ifp);
   *  @endcode
   *
   *  
   *  After running the above code, 
   *  the value of <tt>gI.molarVolumes</tt> is [1.3E-5, 1.3E-5, 2.0E-5].
   * 
   *  Also, we may alternatively query the <tt>besmd</tt> structure to find the
   *  the value in two ways:
   *
   *  @code
   *   LineEntry *le = besmd->searchLineEntry("Molar Volumes");
   *   const double * molarVols = *(static_cast<const double *>(le->currentValueAsVoidP()));
   *   int vecLength = le->get_NumTimesProcessed();
   *
   * 
   *   LineEntry *le = besmd->searchLineEntry("Molar Volumes");
   *   LE_VecDbl *le_dbl = dynamic_cast<LE_VecDblVarLength *>(le);
   *   const double * molarVols = le_dbl->currentTypedValue();
   *  @endcode
   *
   *  Using this alternative, we could have set the external address used in
   *  the constructor of <B>%LE_VecDbl</B> to 0.
   *
   *
   *  <H2> Setting Dependencies on this card </H2>
   *
   *  Dependencies that this card may have on other cards in the input deck
   *  may be specified by calls to the member function #declareDependency().
   *  Dependencies that may be set on this card are listed below:
   *
   *       -  <b>BIDRT_PT_ERROR</b>
   *                 The target %BaseEntry must have occurred before this 
   *                 card in the input deck.
   *       -  <b>BIDRT_ANTITHETICAL_ERROR</b>
   *                 The target %BaseEntry must not occur at all
   *                 in the input deck if this card occurs
   *       -  <b>BIDRT_ZERONUMTIMESREQUIRED</b>
   *                 If the dependency check on the target %BaseEntry
   *                 is satisfied, then the number of required
   *                 times this card is needed is turned to zero.
   *       -  <b>BIDRT_ONENUMTR</b>
   *                 If the dependency check on the target %BaseEntry
   *                 is satisfied, then the number of required
   *                 times this card is needed is turned to one.
   *       -  <b>BIDRT_USETOPROCESS </b>
   *                 The target %BaseEntry must have been already been
   *                 processed, and the default value for this card
   *                 is set to the value of the target %BaseEntry.
   *       -  <b>BIDRT_RTINTMM_ERROR</b>
   *                 If the #CurrValue of this entry is
   *                 in a specified range of values, then a
   *                 dependency check is made against the target
   *                 %BaseEntry.
   *
   *  This card may service dependency requests from other cards
   *  using the base service request, #ansDepCheck().
   *
   * @ingroup blockentryModule
   */
  class LE_VecDblVarLength : public LineEntry {
  public:

    //! Main constructor for a double keyline entry.
    /*!
     *   This sets up the line entry special case.
     *   We make sure to call the base class constructor here to do
     *   much of the initialization. 
     *
     *  When the keyline in the input file of the following form is
     *  found:
     *
     *  "KeyName" = [boolean]
     *
     * The double value at address, addrVal, is assigned the value
     * read in from the input file.
     *
     * @param keyName   C character string setting up the name
     *                  of the keyline to match
     * @param hndVec    Address of the double *, external to the 
     *                  object, which will get assigned the vector
     *                  of doubles
     * @param vecLength initial Length of the vector of doubles.
     *                   (it may change)
     * @param numRL     Number of Required lines in the input file.
     *                  A fault is triggered if this number is nonzero
     *                  and the keyline isn't found in the input file.
     * @param varName   Variable name that is defined by this command.
     *                  This is only used for IO purposes.
     */
    LE_VecDblVarLength(const char *keyName, double **hndVec, int vecLength,
		       int numRL = 0, const char *varName = 0);

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    LE_VecDblVarLength(const LE_VecDblVarLength& right);

    //! Copy assignment operator
    /*!
     * @param right Object to be copied
     */
    LE_VecDblVarLength& operator=(const LE_VecDblVarLength& right);

    //! Duplicator function
    /*!
     * This function duplicates the entry and returns a pointer
     * to a LineEntry
     */
    virtual LineEntry* duplMyselfAsLineEntry() const;

    //! Destructor
    virtual ~LE_VecDblVarLength();

 
    //! Process this line Entry, which assigns a double variable
    /*!
     * This function is called when it has been determined that
     * the current KeyName from the input file matches this 
     * object's KeyName.
     *
     * This function then processes the entry.
     * Processing involves checking for the satisfaction 
     * of runtime dependencies.
     * Then, the CurrentValue is assigned and the external pointer
     * to boolean is assigned the CurrentValue.
     * 
     * The NumTimesProcessed field is also incremented.
     *
     * @param lineArgTok  Pointer to the token containing the
     *                    arguments to the lineEntry. This is 
     *                    everything after the "=" sign.
     */
    void process_LineEntry(const TK_TOKEN *lineArgTok);

    //!  This routine will print out a processed line
    /*!
     *  The default behavior is to print the original line with a "====>"
     *  prefix to indicate that action has been taken on it.
     *
     * @param lineArgTok  Pointer to the token containing the
     *                    arguments to the lineEntry. This is 
     *                    everything after the "=" sign.
     */
    void print_ProcessedLine(const TK_TOKEN *lineArgTok) const;

    //! Print out API information about this keyline
    /*!
     * This routine will print out to stdout information about
     * the keyline. This command is used to document the interface.
     *
     * @param ilvl Level of the indentation to use in printing
     */
    void print_usage(int ilvl) const;


    //! Adjust base address to store the value
    /*!
     * @param addrAdjustment Offset in raw bytes to adjust the
     *              internal value of the Address that will receive
     *              the boolean.
     */
    virtual void adjustAddress(LONG_PTR addrAdjustment );

  //! Return the current vector of values as a const double *
    /*!
     * This is a nonvirtual function since the return type
     * is specific to this child.
     */
    const double * currentTypedValue() const;

    //! Return the current value as a const pointer to void
    /*!
     *  The calling function must know how to interpret the
     *  pointer to void. This is not very far fetched, since
     *  it should know that the value must be a double.
     *
     *  @return Returns a pointer to current vector of values
     *          internally kept by the object.
     */
    virtual const void * currentValueAsVoidP() const;

    //! Set the default value of this card
    /*!
     * @param defValue Boolean value that will be the default.
     */
    void set_default(double defValue);

    //! Check for requirements being met at the end of input
    /*!
     * @param throwSpecificError Usually false.
     */
    virtual bool checkRequirements(bool throwSpecificError = false);

    //!  Set the maximum and minimum value of an entry
    /*!
     *  An error will be thrown inf the keyline value isn't
     *  between these two.
     *
     * @param maxValue maximum value.
     *           The default value is DBL_MAX
     * @param minValue minimum value
     *           The default value is -DBL_MAX
     */
    void set_limits(double maxValue, double minValue);

    //! Sets the string that will be used as the name of 
    //! the variable.
    /*!
     *  @param ps  string name of the variable.
     */
    void set_PrintString(const char *ps);

  private:
    
    //! We store the handle to the double vector which will be
    //! filled with the output.
    double **HndlDblVec;

    //! Max value that this number can attain
    double MaxVal;

    //! Min value that the entry can attain
    double MinVal;

    //! Default value
    double DefaultVal;

    //! Length of the vector to be filled
    int    VecLength;

    //! current index of the internal vector
    int    CurrIndex;

    //! Current value
    /*!
     *  Initially this gets set to the default. Then when
     * the keyline is called, this gets set to the value of 
     * the argument of the keyline.
     */
    double CurrValue;

    //! Vector of entered values.
    /*!
     * Length is the original length of the charList.
     * This is initially set to the default values.
     */  
    double *m_CurrentVecValues;

    //! PrintString -> name of the variable that this keyLine update,
    //! as it will appear on descriptive printouts.
    /*!
     * This defaults to the name of the keyline
     */
    char PrintString[MAX_INPUT_STR_LN+1];
  };
}
#endif
