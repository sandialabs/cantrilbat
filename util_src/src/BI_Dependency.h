/**
 * @file BI_Dependency.h
 *   Declarations for the BI_Dependency base class
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
#ifndef BI_DEPENDENCY_H
#define BI_DEPENDENCY_H

#include "BI_InputError.h"

namespace BEInput {
  /*
   *
   * DEPENDENCIES
   * --------------------
   *
   * A line entry or block entry may have multiple dependencies attached to
   * the entry. These entries belong to the line or block entry. However, they point to
   * a different line or block entry which contains information or is a
   * prerequisites for the current entry to be meaningfull.
   *
   * A dependency has two parts. The first part determines whether the
   * dependency is satisfied. The type of the dependency is specified by
   * the BIDT_TYPE enum value.  The second part determines what action
   * to be taken when the dependency is satisified or when the dependency
   * is not satisfied. The action is set by the BIDRT_TYPE enum value.
   *
   * For example if the initial guess for a simulation is to be obtained 
   * from several methods, then the second card "Saved Guess file Name"
   * only makes sense if the Initial Guess field states File and not Zero.
   *
   * case 1:
   *     Initial Guess = File
   *     Saved Guess file Name = init.txt
   *
   * case 2:
   *     Initial Guess = Zero
   *     # Saved Guess file Name = init.txt
   *
   * case 3: (Error Condition)
   *     Initial Guess = Zero
   *     Saved Guess file Name = init.txt
   *
   * case 4: (Error Condition)
   *     Initial Guess = File
   *     # Saved Guess file Name = init.txt
   *
   * case 5: (Error Condition - wrong order)
   *     Saved Guess file Name = init.txt
   *     Initial Guess = File

   *
   * We would like case 1 and case 2 to be valid inputs. But, we would 
   * like cases 3 and 4 to produce an error condition, to avoid ambiguity in 
   * the mind of the user.
   *
   * In order to accomplish this we want 2 things. First, We require that
   * the "Initial Guess" card be in the input deck and be processed first
   * before the "Saved Guess" card.
   * Second, we require that the "Initial Guess" card have a certain value.
   *
   * Also, we may want to require ordering wihin an otherwise free-form
   * input deck, for the sake of consistency and readibility. Therefore,
   * case 5, we will want to consider an error as well.
   * 
   * Thus, we define a dependency check, attaching it to the second card,
   * such that it will ask the first card, "Initial Guess", what value 
   * it has. The dependency will know what value the Initial Guess card
   * needs in order for the dependency/prerequisite to be satisfied.
   *
   * Below is the code necessary to set up the reading of the input
   * deck with the required error conditions flagged. Note, the
   * userInput struct is set up in global scope so that it stays
   * in existence long enough.
   *
   * @code
   *   struct userInput {
   *      userInput() :
   *        m_initialGuessType(1),
   *        m_fileName("init.txt") {};
   *      int m_initialGuessType;
   *      str::string m_fileName;
   *   };
   *   userInput uI;
   *   BlockEntry *beamd = new BlockEntry();
   *
   *    // Set up the Initial Guess card and make it mandatory
   *   char *cuse[2] = {"file", "zero"};
   *   LE_PickList *lep_ig = new LE_PickList("Initial Guess",
   *                                         &(uI->m_initialGuessType),
   *                                         cuse, 2, 1, "m_initialGuessType");
   *   lep_ig->set_default("zero");
   *   beamd->addLineEntry(lep_ig);
   *    
   *   LE_OneStr *les_gf = new LE_OneStr("Saved Guess File Name",
   *                                     &(uI->m_fileName),
   *                                     1, 1, 0, "SavedGuessFileName");
   *   les_gf->set_default("init.txt");
   *   beamd->addLineEntry(les_gf);
   *
   *   // If Initial Guess is a file, make the "Saved Guess File Name" card mandatory
   *   // We do this with the result type BIDRT_ONENUMTR, which indicates that the
   *   // current card must exist one time in the solution deck.
   *   BI_DepIntMaxMin* depimm_ig =
   *               new BI_DepIntMaxMin(lep_ig, BIDT_INTMAXMIN,
   *                                   0, 0, BIDRT_ONENUMTR);
   *   les_gf->declareDependency(depimm_ig);
   *
   *   // If Initial Guess is a zero make the appearance of this card in the
   *   // input deck an error.
   *   BI_DepIntMaxMin* depimm_ig2 =
   *               new BI_DepIntMaxMin(lep_ig, BIDT_INTMAXMIN,
   *                                   1, 1, BIDRT_ANTITHETICAL_ERROR);
   *   les_gf->declareDependency(depimm_ig2);
   *
   *   // Make sure the Initial Guess card occurs before File Name card
   *   les_gf->declareSimpleOrderDependency(lep_ig);
   *
   * @endcode
   *
   *
   *  Listing of when Dependencies are checked:
   *
   *  BlockEntries:
   *    When a BlockEntry is encountered in the input deck,  it's list of dependencies is 
   *    checked for BIDRT_PT_ERROR ResultTypes. It checks for satisfication
   *    of these dependencies, and error exists if it finds one that
   *    is not satisfied.
   *
   *    At the start or processing a BlockEntry encountered in the input deck
   *    (see BlockEntry::Initialization(), it's list of dependencies is
   *    checked for BIDRT_ANTITHETICAL_ERROR result types.
   *    It checks for satisfaction of these dependencies, and error
   *    exists if it finds one that is satisfied.
   *
   *    After the entire input deck is processed (see BlockEntry::checkRequirements()),
   *    it's list of dependencies is 
   *    checked for BIDRT_ZERONUMTIMESREQUIRED and BIDRT_ONENUMTR ResultTypes. 
   *    It checks for satisfication of these dependencies, and modifies
   *    its own value of n_numTimesRequired based on the output.
   *
   *    After the entire input deck is processed (see BlockEntry::checkRequirements()),
   *    it's list of dependencies is 
   *    checked for BIDRT_ANTITHETICALERROR ResultTypes. 
   *    It checks for satisfication of these dependencies, and
   *    error exits if it finds one.
   *
   *  LineEntry:
   *    At the start or processing a LineEntry encountered in the input deck
   *    (see LineEntry::process_LineEntry(), it's list of dependencies is
   *    checked for BIDRT_PT_ERROR result types.
   *    It checks for satisfication of these dependencies, and error
   *    exists if it finds one that is not satisfied.
   *
   *    At the start of processing a LineEntry encountered in the input deck,
   *    (see LineEntry::process_LineEntry(), the list of dependencies for that LineEntry is
   *    checked for dependencies with the BIDRT_ANTITHETICAL_ERROR result types.
   *    It checks for satisfication of these dependencies, and error
   *    exists if it finds one that is satisfied.
   *
   *    After the entire input deck is processed (see LineEntry::checkRequirements()),
   *    a LineEntry or BlockEntry's list of dependencies is 
   *    checked for BIDRT_ZERONUMTIMESREQUIRED and BIDRT_ONENUMTR ResultTypes. 
   *    It checks for satisfication of these dependencies, and modifies
   *    its own value of n_numTimesRequired based on the output.
   *
   *    After the entire input deck is processed  (see LineEntry::checkRequirements()),
   *    it's list of dependencies is 
   *    checked for BIDRT_ANTITHETICALERROR ResultTypes. 
   *    It checks for satisfication of these dependencies, and
   *    error exits if it finds one.
   * 
   *  Specific LineEntries:
   *
   *   LE_OneInt, LE_PickList:
   *     At the start of processing a LE_PickList encountered in the input deck
   *    (see LE_PickList::process_LineEntry(), it's list of dependencies is
   *    checked for BIDRT_USEDTOPROCESS result types.
   *    It checks for satisfication of this dependencies, and error
   *    exists if it finds one that is not satisfied. It uses the return
   *    int to set the default value for the LE_PickList.
   *
   *   LE_VecDbl, LE_VecDblVarLength:
   *    At the start of processing a LE_VecDlb encountered in the input deck
   *    (see LE_VecDbl::process_LineEntry(), it's list of dependencies is
   *    checked for BIDRT_USEDTOPROCESS result types.
   *    It checks for satisfication of this dependencies, and error
   *    exists if it finds one that is not satisfied. It uses the return
   *    int to set the default vector length of the vector of doubles
   *    to be read in.
   *
   *   LE_OneStr
   *    At the start of processing a LE_OneStr encountered in the input deck
   *    (see LE_OneStr::process_LineEntry(), it's list of dependencies is
   *    checked for BIDRT_USEDTOPROCESS result types.
   *    It checks for satisfication of this dependencies, and error
   *    exists if it finds one that is not satisfied. It uses the return
   *    int to set the default number of tokens to be expected in the
   *    string to be read in.
   *
   *
   *  LE_PickList, LE_OneBoolInt:
   *    At the end of processing a LE_PickList encountered in the input deck
   *    (see LE_PickList::process_LineEntry(), it's list of dependencies is
   *    checked for BIDRT_RTINTMM_ERROR result types.
   *    It checks for satisfaction of this dependencies, and error
   *    exists if it isn't satisfied. To check for satisfaction of the
   *    dependency, it checks to see whether the current int value
   *    from the LE_PickList process is between a max and min value
   *    set within the dependency.
   *
   *
   *
   */


  //!  enum for specification of the dependency type information 
  /*!
   *  The dependency type is a strict definition of what it takes
   *  for the function checkDependency() to return a true or false condition
   */
  enum BIDT_TYPE {

    //! Target BaseEntry has already been processed
    /*!
     * BIDT_ENTRYPROCESSED -> This dependency states that for an entry to 
     *                        be valid  another specified entry must have 
     *                        appeared before this entry in the input deck. 
     *                        The entry must have appeared even if it is 
     *                        officially an optional entry.
     */
    BIDT_ENTRYPROCESSED = 0,

    //! Target BaseEntry has already been processsed and has/will supply
    //! a single int to the dependent BaseEntry.
    /*!
     * BIDT_ONEINT -> This dependency involves all of what  BIDT_ENTRYPROCESSED
     *                does. And, it requires that the required entry provide
     *                a single integer to the dependent entry at the time that
     *                dependent entry is called to process a line.
     *                Thus, BIDT_ONEINT only makes sense if the required 
     *                BaseEntry
     *                somehow can suppy an int value, and if the dependent BaseEntry
     *                can somehow use that int value.
     */
    BIDT_ONEINT,

    //!  Target BaseEntry has already been processsed and the supplied
    //!  single int is within a max/min int value.
    /**
     * BIDT_INTMAXMIN -> This dependency involves all of what  BIDT_ENTRYPROCESSED
     *                does. And it requires that the required BaseEntry somehow
     *                supply an int value by hook or crook. That int value is
     *                checked to see if it is between a max and min bounds.
     *                If it is, then the dependency is satisfied. If it snot
     *                then the dependency is not satisfied, and an error condition
     *                is thrown.
     */
    BIDT_INTMAXMIN
  };


  //! enum ( BIDRT_#####) specifying the result type to be taken
  //! upon determinination that the dependency check has been satisfied
  //! It also specifies the timing of the dependency check.
  /*!
   *  The default case here is to see if the dependency is satisfied.
   *  If it is, nothing occurs. If it isn't, a BI_InputError is thrown.
   *  However, there are other cases as well.
   */
  enum BIDRT_TYPE {

    //! Default behavior is to throw an error if a dependency isn't satisfied
    /*!
     * The default result is that if a dependency isn't met at the time of
     * processing, an BI_InputError error condition is created and thrown.
     * If the dependency is met, nothing further is done.
     *  (default in most cases)
     * This is used to ensure ordering within the input deck.
     *
     *  timing -> at the the time the dependent BaseEntry is processed, the
     *            dependency must be satisfied.
     */
    BIDRT_PT_ERROR = 0,

    //!  A satisfied dependency results in the dependent BaseEntry becoming optional
    /*!
     *    This result indicates that if the dependency check is met, then
     *    the number of times required for the current BaseEntry is set
     *    to zero. Therefore, it becomes an optional card.
     *
     *  timing -> at the end of processing of the enclosed block,
     *            using checkRequirements()
     */
    BIDRT_ZERONUMTIMESREQUIRED,

    //!  A satisfied dependency results in the dependent BaseEntry 
    //!  becoming mandatory
    /*!
     *    This result indicates that if the dependency check is met, then
     *    the number of times the line entry must be processed is equal
     *    to one. Basically, we have set an optional card to
     *    a required card by this action.
     *
     *  timing -> at the end of processing of the enclosed block,
     *            using checkRequirements()
     */
    BIDRT_ONENUMTR,

    //! A satisfied dependency is needed because the dependent
    //! BaseEntry needs the information to process its own keyline
    //! Information
    /*!
     *    This result indicates that if the dependency check is met, then
     *    the information exchanged is needed for the processing
     *    of the owning line entry at the time the  owning keyline is encountered in
     *    the input deck.
     *    An error condition is set if the dependency is set but not
     *    met when the owning keyline is read.
     *    This is because a previous bit of information is needed
     *    to process the current entry.
     *    This is usually used to set the default value of the owning
     *    keyline's entry (e.g., the picklist value),
     *    based on a previous entry. Or, it's used to set the malloc
     *    length of a vector that needs to be allocated before it is used
     *    (e.g., LE_VecDbl())
     *
     *  timing -> at the the time the dependent BaseEntry is processed.
     */
    BIDRT_USETOPROCESS,

    //! The dependency check is only made if the dependent BaseEntry
    //! contains an int value between a max int and a min int. Then,
    //! a regular dependency check on the target BaseEntry is made.
    /*!
     *     This result indicates that a dependency check of the type
     *     specified by BIDDT_#### will be carried out iff the local
     *     owning BaseEntry processing
     *     produces an integer value between a MAX and MIN value
     *     stored within the dependency check itself.
     *     If so, then a regular dependency check is carried out
     *     on the required BaseEntry.
     *     Thus, for example if a particular PickList entry needs a
     *     certain precursor card to have occurred and to have a
     *     particular value, we can flag that condition with this
     *     dependency condition. Failures will throw an error condition.
     *
     *  timing -> at the the time the dependent BaseEntry is processed.
     */
    BIDRT_RTINTMM_ERROR,

    //! A satisfied dependency causes an immediate error exit
    /*!
     *     This result indicates that if a dependency check
     *     is met on the required item, then an immediate error is
     *     thrown from the dependent item.
     *     This may be used when 2 cards or blocks are mutually
     *     exclusive, ie., if you had that card, then you can't use this card.
     *
     *  timing -> at the the time the dependent BaseEntry is processed.
     */
    BIDRT_ANTITHETICAL_ERROR
  };

  //! Dependency service request type
  /*!
   *       Right now this is figured out from the previous two ints.
   *       It is not part of the interface.
   *
   *     The service request is a format of the request for dependency 
   *     satisfaction that is made between the dependency object and the
   *     required BaseEntry. The required BaseEntry must be able to service
   *     that request. 
   *
   */
  enum BIDSR_TYPE {

    //! Simple Ordering and query about whether it has been called previously
    /*!
     *       This service request just asks the required BaseEntry whether it has
     *       been called previously.
     *           It calls the boolean function requiredEntry->ansDepCheck();
     */
    BIDSR_SIMPLEORDERING = 0,

    //!       This service request asks the required BaseEntry whether it has
    //!       been called previously and asks it to return an int.
    /*!
     *     It calls the boolean function,
     *                  requiredEntry->checkDepOneInt(&int retnInt);
     *       If the required BaseEntry hasn't been called yet, it returns
     *       the default retnInt value. If it has, then it returns the retnInt
     *       specified by the BaseEntry. For example, if the BaseEntry is a 
     *       LE_PickList object, then the integer position of the picked entry
     *       is returned as the integer.
     */
    BIDSR_ONEINT 
  };

  /* forward declaration */
  class BaseEntry;

  //!  This base class handles dependencies between objects.
  //!  This class describes what the dependency for the owning object is.
  /*!
   * A line entry or block entry may have multiple dependencies attached to
   * the entry.
   * These entries belong to a line or block entry. However, they point to
   * a target line or block entry which contains information or is a 
   * prerequisite for the current entry to be meaningfull.
   *
   */
  class BI_Dependency {
  public:

    //! Default Constructor
    /*!
     * @param be         const pointer to the target Base entry that 
     *                   is the source of the dependency
     *                   This must be a valid BaseEntry object
     * @param BIDT_type  Type of dependency
     * @param BIDRT_type Result of the dependency 
     *                   The default is to ask for an ordering 
     *                   dependency that is satisfied at the
     *                   time that the Dependent BaseEntry is
     *                   processed.
     */
    BI_Dependency(const BaseEntry *be, BIDT_TYPE BIDT_type, 
		  BIDRT_TYPE BIDRT_type = BIDRT_PT_ERROR);

    //! Copy constructor
    /*!
     * @param b object to be copied
     */
    BI_Dependency(const BI_Dependency& b);

    //! Assignment operator
    /*!
     *  @param b object to be copied
     */
    BI_Dependency& operator=(const BI_Dependency&b);

    //! Virtual destructor
    virtual ~BI_Dependency();

    //! Base Class Duplicator function
    virtual BI_Dependency* duplicateMyself() const;

    //! Checks to see if the dependency is satisfied
    /*!
     * This function calls ansDepCheck() on the target
     * BaseEntry to determine whether the dependency has
     * been satisfied.
     *
     * @return returns a boolean 
     */
    virtual bool checkDependencySatisfied() const;

    //! Checks to see if the dependency is satisfied
    /*!
     * This function calls ansDepCheckOneInt() on the target
     * BaseEntry to determine whether the dependency has
     * been satisfied. 
     * It stores the returned int m_RetnOneInt
     *
     * @param returnInt Output variable containing the int
     *                  value returned from the target
     *                  BaseEntry
     *
     * @return returns a boolean indicating whether
     *                 the dependency is satisfied
     */
    virtual bool checkDepOneInt(int &returnInt) const;

    //! Check an int against the permissible values
    /*!
     *  Checks the int argument against permissible max
     *  and min values stored in this object
     *
     * @param cardInt input int argument
     *
     * @return returns boolean indicating whether the
     *         argument int is in range.
     */
    virtual bool checkCardIntMaxMin(int cardInt);

    //! Set the max and min int range values
    /*!
     * @param cMax Maximum int value
     * @param cMin minimum int value
     */
    virtual void setCardIntMaxMin(int cMax, int cMin);

    //! Returns the ResultType of this dependency
    BIDRT_TYPE ResultType() const;

    //! Returns the ServiceRequestType of this dependency
    BIDSR_TYPE ServiceRequestType() const;

    //! Returns a const pointer to the target BaseEntry
    const BaseEntry *TargetBaseEntry() const;

  protected:

    //! Pointer to the target BaseEntry for this dependency
    /*!
     * Must be nonnull
     */
    const BaseEntry *RequiredEntry;
  
    //! This is the exact definition of what it means for the
    //!  dependency to be satisfied.
    /*!
     *  The dependency type is a strict definition of what it takes
     *  for the function checkDependency()
     *           BIDT_ENTRYPROCESSED
     *           BIDT_ONEINT 
     *           BIDT_INTMAXMIN 
     */ 
    BIDT_TYPE m_DT_type;

    //! This is the result -> Basically what happens if the
    //! dependency is met.
    /*!
     * The "dependency result" type defines the action to be taken
     * when the dependency is either satisfied or not satisfied.
     *
     *           BIDRT_PT_ERROR
     *           BIDRT_ZERONUMTIMESREQUIRED 
     *           BIDRT_ONENUMTR
     *           BIDRT_USETOPROCESS
     *           BIDRT_RTINTMM_ERROR
     */
    BIDRT_TYPE m_RT_type;

    //!  Dependency service request type
    /*!
     *      -> right now this is figured out from the previous two ints.
     *         It is not part of the interface.
     *          
     */
    BIDSR_TYPE m_SR_type;

    //! Max int value for the Dependent BaseEntry
    /*!
     * These max and min values are used in the BIDRT_RTINTMM_ERROR case.
     * When the local owning BaseEntry entered value is between these values,
     * then a dependency check is made by calling the required BaseEntry object.
     */
    int CardMax;

    //! Min int value for the Dependent BaseEntry
    /*!
     * These max and min values are used in the BIDRT_RTINTMM_ERROR case.
     * When the local owning BaseEntry entered value is between these values,
     * then a dependency check is made by calling the required BaseEntry object.
     */
    int CardMin;
  };
}
#endif
