/**
 * @file LE_OneDblUnits.h
 *  Header for the LineEntry of a single double with units conversion
 *  (see \ref blockentryModule and class 
 *  \link BEInput::LE_OneDblUnits LE_OneDblUnits\endlink).
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef LE_ONEDBLUNITS_H
#define LE_ONEDBLUNITS_H

#include "LE_LineEntry.h"

namespace BEInput {
  class BE_UnitConversion;

  /*!
   *   This sets up the <b>LineEntry</b> special case for the entry of a single
   *   double input with units conversion.  It declares interfaces for specifying the
   *   entry and for processing the entry once a match is
   *   made in the input deck. The general form of the input is
   *
   *   \f[
   *     KeyName = \mathrm{double\qquad} \left[ \mathrm{units\_string} \right]
   *   \f]
   *
   *  The <I>KeyName</I> variable  may be consist of multiple
   *  tokens. It is tokenized and compared in a case insensitive format
   *  at the start of the processing and before any matching
   *  is carried out.
   *  \f$ \mathrm{double} \f$ is a single double. The code will check that
   *  the rhs consists of one and only one well-formed double.
   *  Additionally, the following tokens are understood and accepted in
   *  place of \f$ \mathrm{double} \f$ : <tt>DBL_MAX</tt>,  <tt>DBL_MIN</tt>,
   *  and <tt>Default</tt>. If <tt>DBL_MAX</tt> or  <tt>DBL_MIN</tt> is
   *  encountered, the value is set to the respective value as defined
   *  in the file math.h. If the value <tt>Default</tt> is encountered,
   *  the #CurrValue is set to the default value, #DefaultVal,
   *  specified for the object.
   *
   *  The \f$ \mathrm{units\_string} \f$ is an optional string. If present
   *  it will trigger a units conversion process to an internally
   *  storred MKS units value. Note, the dimensionality of the units string
   *  is checked via specification of the expected dimensionality
   *  at the time the line input is set up. For example, if a single
   *  length is expected, then the \f$ \mathrm{units\_string} \f$
   *  value must correspond to a string representing a single length.
   *  See BEInput::BE_UnitConversion for more information on this process.
   *
   *  If the default value specified for the object 
   *  is set to <tt>NO_DEFAULT_DBL</tt>, which is the initial value for
   *  #DefaultVal, then an error is thrown when <tt>default</tt> 
   *  is encountered on the rhs.
   *
   *  When a match is found, the rhs of the card is processed into a double
   *  using the function #process_LineEntry().
   *  The double is then assigned to the #CurrValue of this object.
   *  And, it is optionally written out to an external address of an 
   *  int, #AddrVal,
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
   *   the global address globInput.viscosity.
   *
   * @code
   *      struct globInput {
   *         double Pressure;
   *      } gI;
   *
   *      BlockEntry *besmd = new BlockEntry("Fluid Properties");
   *      int reqd = 1; 
   *      BE_UnitConversion *ucPres = BE_UnitConversionPressure::units();
   *      LE_OneDblUnits *d2 = new LE_OneDblUnits("Pressure",
   *                                              &gI.Pressure,
   *                                              reqd, "Pressure", ucPres);
   *      d2->set_default(101325.);
   *      d2->set_limits(10000000., 1.0E-4);
   *      besmd->addLineEntry(d2);
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
   *         Pressure = 1.0 atm
   *       End Block   Fluid Properties
   *  @endcode
   *
   *  It is an error for the Block "Fluid Properties" to not
   *  have exactly one "Pressure" %LineEntry.
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
   *  the value of <tt>gI.Pressure</tt> is 101325.
   * 
   *  Also, we may alternatively query the <tt>besmd</tt> structure to find the
   *  the value in two ways:
   *
   *  @code
   *   LineEntry *le = besmd->searchLineEntry("Pressure");
   *   double value = *(static_cast<const double *>(le->currentValueAsVoidP()));
   * 
   *   LineEntry *le = besmd->searchLineEntry("Pressure");
   *   LE_OneDbl *le_dbl = dynamic_cast<LE_OneDblUnits *>(le);
   *   double value = le_dbl->currentTypedValue();
   *  @endcode
   *
   *  Using this alternative, we could have set the external address used in
   *  the constructor of <B>%LE_OneDblUnits</B> to 0.
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
   *  using the base service request, #ansDepCheck(), and the service request that
   *  returns an int value, #ansDepCheckOneInt(). The int value returned is the
   *  converted int value of #CurrValue of the %LE_OneDbl object.
   *
   * @ingroup blockentryModule
   */
  class LE_OneDblUnits : public LineEntry {
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
     * @param addrVal   Address of the double, external to the 
     *                  object, which will get assigned the value of 
     *                  the expression. (default 0)
     * @param numRL     Number of Required lines in the input file.
     *                  A fault is triggered if this number is nonzero
     *                  and the keyline isn't found in the input file.
     * @param varName   Variable name that is defined by this command.
     *                  This is only used for IO purposes.
     *
     * @param uc        Valid units conversion object that will be
     *                  used to convert the input double. This object
     *                  is owned by the LE_OneDblUnits object, and will
     *                  deleted during the destructor.
     */
    LE_OneDblUnits(const char *keyName, double *addrVal, int numRL = 0,
		   const char *varName = 0, BE_UnitConversion *uc = 0);

   //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    LE_OneDblUnits(const LE_OneDblUnits& right);

    //! Copy assignment operator
    /*!
     * @param right Object to be copied
     */
    LE_OneDblUnits& operator=(const LE_OneDblUnits &right);

    //! Duplicator function
    /*!
     * This function duplicates the entry and returns a pointer
     * to a LineEntry
     */
    virtual LineEntry* duplMyselfAsLineEntry() const;

    //! Destructor
    virtual ~LE_OneDblUnits();

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
    virtual void adjustAddress(LONG_PTR addrAdjustment);

  //! Return the current value as a double
    /*!
     * This is a nonvirtual function since the return type
     * is specific to this child.
     */
    double currentTypedValue() const;

    //! Return the current value as a const pointer to void
    /*!
     *  The calling function must know how to interpret the
     *  pointer to void. This is not very far fetched, since
     *  it should know that the value must be a double.
     *
     *  @return Returns a pointer to current double value as a
     *          pointer to void.  
     */
    virtual const void * currentValueAsVoidP() const;
  
    //! Set the default value of this card
    /*!
     * @param defValue Boolean value that will be the default.
     */
    void set_default(double defValue);

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

    //!The address of the double that gets updated by this command
    /*!
     *  If this is zero, then no external update is carried out.
     */
    double *AddrVal;

    //! Max value that this number can attain
    double MaxVal;

    //! Min value that the entry can attain
    double MinVal;

    //! Default value
    double DefaultVal;

    //! Current value
    /*!
     *  Initially this gets set to the default. Then when
     * the keyline is called, this gets set to the value of 
     * the argument of the keyline.
     */
    double CurrValue;

    //! PrintString -> name of the variable that this keyLine update,
    //! as it will appear on descriptive printouts.
    /*!
     * This defaults to the name of the keyline
     */
    char PrintString[MAX_INPUT_STR_LN+1];

    //! Pointer to the units conversion object
    BE_UnitConversion *m_uc;
  };
}

#endif
