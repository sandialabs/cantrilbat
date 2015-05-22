/**
 * @file BE_BlockEntry.h
 *    Declarations for BlockEntry object
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef BLOCKENTRY_H
#define BLOCKENTRY_H

/*
 * Need to get the definition of the token structure
 */
#include "BaseEntry.h"

#include <set>

namespace BEInput {

  // Forward declaration
  class LineEntry;

#define BE_ANY_INDEX -7777779

  //!  The BlockEntry Class is used to describe blocks 
  /*!
   *   A block has the following structure
   *
   *    Start Block KeylineName = "start arguments"
   *       . . .
   *       keylines and blocks
   *       . . .
   *    End Block KeylineName = "end arguments"
   *
   *  The keylines for the start and end blocks must match
   *  up. It is an error if they don't.
   *
   *  It may be set as being required in the input block or as being optional
   *  in the input block. Keylines within the block may be set as being
   *  required even if the block is optional.
   *
   *  The arguments to the block appear after an equals sign. They are tokenized
   *  and then available for later processing within the block.
   *
   *  This class contains the methods for recursively scanning
   *  nested blocks and reading their input.
   *
   * @ingroup blockentryModule
   */
  class BlockEntry : public BaseEntry {

  public:

    //! Constructor
    /*!
     *   All parameters but the block name default values. The Block name
     *   must be specified
     *
     * @param blockName         Keyline name of the block
     * @param numTimesRequired  Number of times this block is required
     *                          to be found in the input file
     * @param ParentBlock_input   Pointer to the parent of the current block
     *                          default is 0, indicating that this block is the
     *                          top main block that has no parent blocks.
     */
    explicit BlockEntry(const char *blockName,
			int numTimesRequired = 0,
			BlockEntry *ParentBlock_input = 0);

    //! Copy constructor
    /*!
     * @param right object to be copied
     */
    BlockEntry(const BlockEntry& right);

    //! Assignment operator
    /*!
     *  @param right Object to be copied
     */
    BlockEntry& operator=(const BlockEntry& right);

    //! duplicator function
    /*!
     * Duplicates the object and all underlying objects
     * and returns a pointer to BaseEntry
     */
    virtual BaseEntry* duplMyselfAsBaseEntry() const;

    //! duplicator function
    /*!
     * Duplicates the object and all underlying objects
     * and returns a pointer to BlockEntry
     */
    virtual BlockEntry* duplMyselfAsBlockEntry() const;

    //! Destructor
    virtual ~BlockEntry(void);

    //! Remove all line entries and subblocks
    void clear();
  
    //!  read_block handes the I/O of block commands. It is designed so that
    //!  it can be called recursively.
    /*!
     *  The basic idea is that a line of input is read, stripping comments
     *  so that a line must be interpretted.
     *  Then, the line is either interpretted as a keyline or as the
     *  start of another block.
     *
     *  Block starts are indicated by the first two tokens on a 
     *  line being "START BLOCK". Then, the name of the block is
     *  next. Then, an optional "=" sign and an argument list to the
     *  block is next.
     *
     *  Blocks end with lines indicated by the two tokens "END BLOCK",
     *  followed by the name of the block, followed by
     *  an optional "=" sign and an argument list.
     *
     *  Within a block there may be keylines. A keyline is a string
     *  of arbitrary size followed by the "=" key, or a newline, 
     *  whichever comes first. After the "=" sign is the argument to
     *  the keyline.
     *    
     *    Input
     *  --------------
     * @param input_file File pointer to read additional keylines
     *                   in the recursive calls to read_block
     *                   Default = stdin
     *
     * @param startArgPtr Ptr to the token representing the argument
     *                    to the start block call that initiated the block
     *
     * @param parentBlock Pointer to the block class that is the parent
     *                    to the current block. (default of no parent indicates
     *                    that this is the top block)
     *
     *    Output
     *   ------------
     * @param endArgPtr   Pointer to the TOKEN that represents the 
     *                    arguments to the end block line for the current
     *                    block. This is returned to the calling block
     *                    and represents a way to communicate with the
     *                    calling block
     *                    Default is the null address. 
     *
     */
    void read_block(FILE *input_file, TK_TOKEN *endArgPtr, 
		    const TK_TOKEN *startArgPtr, BlockEntry *parentBlock);

    //!  This function will skip ahead to read a subblock within a block
    //!  structure.
    /*!
     *  This function will skip ahead until it matches the START BLOCK
     *  keyline for the current block. It will then start to process
     *  the current block regularly using the read_block() function.
     *
     * @param input_file File pointer to read keylines
     *                   in the recursive calls to read_block().
     *                   It should be positioned above the the START
     *                   BLOCK keyline for the current block.
     *                   Default = stdin
     */
    void read_this_block_only(FILE *input_file);

    //!  skip_block skips a block and all embedded blocks in the input file
    /*!
     *  On entry the file position should be just after the 
     *  block header of the block to be skipped.
     *
     *    Input
     *  --------------
     * @param input_file File pointer to read additional keylines
     *                   in the recursive calls to read_block
     *                   Default = stdin
     *
     * @param blockName  ptr to the keyline token representing the name
     *                   of the block to skip.
     *
     * @param parentBlock pointer to the block class that is the parent
     *                    to the current block. (default no parent indicates
     *                    that this is the top block.
     */
    void skip_block(FILE *input_file, TK_TOKEN *blockName,
		    BlockEntry *parentBlock);

    //! Virtual function called at the start of internally processing
    //! the block
    /*!
     *  This function may be used to start the process of setting up
     *  internal data functions when the block is first called. 
     *  This is also where the current block processes the
     *  arguments specified on the START BLOCK line. 
     *
     *  The default behavior is listed below.
     *
     *   -  An Error exit will occur if 
     *      a blockARgTok string is actually supplied.
     *
     *   -  We increment the NumProcessedBlocks counter, which
     *      keeps track of the number of times this particular object has
     *      been called for parent block.
     *
     *   -  A dependency check is made to make sure that the 
     *      required dependencies for this block have been satisfied.
     *
     * Derived classes may override the default behavior. Usually
     * derived classes will call the base class Initialization function
     * and perhaps do some other processing.
     *
     *  @param ifp_input  File pointer to read additional keylines
     *                    in the recursive calls to read_block
     *                    Default = stdin
     *
     *  @param blockArgPtr pointer to the TOKEN structure representing
     *                     the argument list for the START BLOCK
     */
    virtual void Initialization(FILE *ifp_input, 
				const TK_TOKEN *blockArgPtr);

    //! This virtual function is used to wrap up the setup of the current
    //! block before returning to the parent block.
    /*!
     *
     *  The default behavior is listed below.
     *
     *   -  An Error exit will occur if 
     *      a blockArgTok string is actually supplied.
     *
     *   -  A checkRequirements check is made on all enclosed line
     *      entry objects. If they haven't been satisfied, an error
     *      is generated. 
     *
     *   -  A checkRequirements check is made on all enclosed blocks.
     *      If they haven't been satisfied, an error is generated. 
     *
     * Derived classes may override the default behavior. Usually
     * derived classes will call the base class Initialization function
     * and perhaps do some other processing.
     *
     *  @param ifp_input  File pointer to read additional keylines
     *                    in the recursive calls to read_block
     *                    Default = stdin
     *
     *  @param blockArgPtr pointer to the TOKEN structure representing
     *                     the argument list for the START BLOCK
     */
    virtual void Wrapup(FILE *ifp_input, const TK_TOKEN *blockArgPtr);

    //! Adjust the address of external target objects
    /*!
     *  External addresses are adjusted by a raw byte value.
     *  Note, the input is a incremental adjustment. The cumulative
     *  adjustment is kept within the BlockEntry class object.
     *
     * @param adjustAAA  Raw bytes to adjust all target addressses
     *                   for objects within the block.
     *
     * @note this puts a severe restriction on the type of
     *       object to be used for the target.
     */
    virtual void adjustAddress(LONG_PTR adjustAAA);

    //! Set the default value of this card
    /*!
     * This will trip an error condition, unless overriden in
     * a derived class
     *
     * @param defValue Boolean value that will be the default.
     */
    virtual void set_default(double defValue);

    //! This call zeroes the numTimesProcessed count
    /*!
     * This is needed 
     * when BlockEntry is used in a multiple pass format.
     * It zeroes the numTimesProcessed field for the
     * current block and all enclosed blocks and LineEntry's.
     */
    void ZeroLineCount();

    //! This routine will print a summary of the permissible line Entries
    //! and subblocks for this particular block
    /*!
     *  This is called recursively to write out the entire structure.
     *  It writes to stdout.
     *
     *  @param indent_lvl Current value of the indent lvl.
     */
    virtual void print_usage(int indent_lvl = 0) const;

    //!   Match the block name as defined by a sequence of tokens.
    //!    Return the address of the block entry.
    /*!
     *    This function will lop off unwanted trailing "block" "end" and "start"
     *    words from the token list.
     *
     * @param keyLinePtr Pointer to the TOKEN to be used to match
     *                   the block.
     * @param contribIndex If the contribIndex is 0 or pos, an exact match with the index is
     *                     required. if contribIndex is < 0, the following procedure is used.
     *                     The value of the lowest  multiContribIndex value for a block that
     *                     hasn't been processed gets selected.
     *  @param includedMatch Bookean If true than the TOKEN must only be included in the
     *                       Block name. If false, then an exact match is required
     */
    BlockEntry *match_block(const TK_TOKEN *keyLinePtr,
			    int contribIndex = BE_ANY_INDEX, 
			    bool includedMatch = false) const;

    //!   Match the block name as defined by a sequence of tokens and optionally
    //!   defined by the argument name of the block.
    //!    Return the address of the block entry.
    /*!
     *    This function will lop off unwanted trailing "block" "end" and "start"
     *    words from the token list.
     *
     * @param keyLinePtr Pointer to the TOKEN to be used to match
     *                   the block.
     * @param contribIndex If the contribIndex is 0 or pos, an exact match with the index is
     *                     required. if contribIndex is < 0, the following procedure is used.
     *                     The value of the lowest  multiContribIndex value for a block that
     *                     hasn't been processed gets selected.
     *  @param includedMatch Bookean If true than the TOKEN must only be included in the
     *                       Block name. If false, then an exact match is required
     *  @param keyArgNamePtr If nonnull, then the name is matched as well. 
     */
    BlockEntry *match_block_argName(const TK_TOKEN *keyLinePtr,
                            bool includedMatch = false,
                            int contribIndex = BE_ANY_INDEX,
                            const TK_TOKEN *keyArgNamePtr = 0) const;


    //!   Match the block name as defined by a sequence of tokens.
    //!    Return the address of the block entry.
    /*!
     *    This function will lop off unwanted trailing "block" "end" and "start"
     *    words from the token list.
     *
     * @param keyLinePtr Pointer to the character string to be used to match
     *                   the block.
     * @param contribIndex If the contribIndex is 0 or pos, an exact match with the index is
     *                     required. if contribIndex is < 0, the following procedure is used.
     *                     The value of the lowest  multiContribIndex value for a block that
     *                     hasn't been processed gets selected.
     *  @param includedMatch Bookean If true than the TOKEN must only be included in the
     *                       Block name. If false, then an exact match is required
     */
    BlockEntry *match_block(const char *keyLinePtr,
			    int contribIndex = BE_ANY_INDEX, bool includedMatch = false) const;
  
    //! Do initialization of the subblock at the current block lvl.
    /*!     
     *  Input
     * ------
     *  @param  subBlockPtr = pointer to the subblock object 
     *  @param  keyArgTok   = arguments to the subblock call.
     *
     *  Note, this functionality is currently unused.
     */
    virtual void Initialize_SubBlock(BlockEntry *subBlockPtr,
				     const TK_TOKEN *keyArgTok);


    //! Do wrapup work for the subblock at the current block lvl.
    /*!     
     * 
     * @param subBlockPtr     = pointer to the subblock object 
     * @param subBlockEndPtr  = arguments to the subblock call.
     */
    virtual void Wrapup_SubBlock(BlockEntry *subBlockPtr,
				 const TK_TOKEN *subBlockEndPtr);
 
    //!   Match a LineEntry in the current block
    /*!
     *    See if a match exists between the character string and the keyLine
     *    structure. If it does match, then pick the best match via matching
     *    via the contribIndex. keylines, contribIndex pairs are unique.
     *    For unspecified contribIndex, return the highest value multiContribIndex
     *    value item that matches the lineEntry.
     *
     *    This function doesn't search recursively through subblocks.
     *
     *    @param lineEntryTok  TOKEN ptr for the keyline to be matched.
     *    @param contribIndex  contribIndex value. Defaults to  BE_ANY_INDEX,
     *                         which doesn't match any index.
     *
     *    @return Return a pointer to the LineEntry object that matches.
     */
    LineEntry *match_keyLine(const TK_TOKEN *lineEntryTok,
			     int contribIndex = BE_ANY_INDEX) const;

    //! After a keyline match, the keyline processes the argument list
    /*!
     *  The basic approach is to pass the control of dealing with line entries
     *  down to the individual LineEntry objects themselves.
     *
     *  @param curLE        Pointer to teh current line entry 
     *  @param keyArgTokPtr TOKEN ptr for the argument list
     */
    virtual void process_LineEntry(LineEntry *curLE,
				   const TK_TOKEN *keyArgTokPtr);

    //! Function is called from the owning BlockEntry 
    //! on every LineEntry that doesn't match an object
    /*!
     *  Currently, this is blank
     *
     *  @param ifp    File pointer
     *  @param keyLineTok TOKEN ptr for keyline that doesn't match
     *  @param keyArgTokPtr TOKEN ptr for arguments to the keyline
     */
    virtual void skip_lineEntry(FILE *ifp, 
				const TK_TOKEN *keyLineTok,
				const TK_TOKEN *keyArgTokPtr) const;
  
    //!   This adds a subblock entry structure to the current BlockEntry
    //!   structure.
    /*!
     * @param sb pointer to the subblock BlockEntry object 
     */
    void addSubBlock(BlockEntry *sb);

    //!   This adds a subblock entry structure to the current BlockEntry
    //!   structure.
    /*!
     * @param le pointer to the LineEntry object 
     */
    void addLineEntry(LineEntry *le);
 
    //! Report the number of times a line entry has been processed.
    /*!
     *   This does a recursive search for a Line Entry on the current block
     *   and all subblocks of the current block, using a character string
     *   as the input.
     *
     *  @param lineName Character string reprsentation of the keyLine
     */
    int reportNumProcessedLines(const char *lineName) const;

    //! Check for all requirements being met at the end of input
    //! for that block
    /*!
     * @param throwSpecificError If true then you should throw a 
     *        specific error condition, if you have one. If not,
     *        an generic error condition will be thrown on return.
     */
    virtual bool checkRequirements(bool throwSpecificError);


    //!  This does a recursive search for a Line Entry on the current block
    //! and all subblocks of the current block
    /*!
     * It uses a character string as the input.
     *
     * @param lineName character string as input
     */
    LineEntry *searchLineEntry(const char * const lineName) const; 

    //!  This does a recursive search for a Line Entry on the current block
    //! and all subblocks of the current block
    /*!
     * It uses a TOKEN ptr as the search input
     *
     * @param nameLE character string as input
     */
    LineEntry *searchLineEntry(const TK_TOKEN * const nameLE) const;

    //! This does a recursive search for a Block Entry 
    //! name under the current block
    //!  and under all subblocks of the current block.
    /*!
     * It uses a character string as the input.
     *
     * @param bName block name
     * @param includedMatch Bookean If true than the TOKEN must only be included in the
     *                       Block name. If false, then an exact match is required
     */
    BlockEntry *searchBlockEntry(const char * const bName, bool includedMatch = false,
                                 int contribIndex = BE_ANY_INDEX, const TK_TOKEN * const blockArgName = 0) const;

    //! This does a recursive search for a Block Entry 
    //! name under the current block
    //! and under all subblocks of the current block.
    /*!
     * It uses a TOKEN ptr as the search input
     *
     * @param nameBN block name
     * @param includedMatch Bookean If true than the TOKEN must only be included in the
     *                       Block name. If false, then an exact match is required
     */
    BlockEntry *searchBlockEntry(const TK_TOKEN * const nameBN, bool includedMatch = false, 
                                 int contribIndex = BE_ANY_INDEX, const TK_TOKEN * const blockArgName = 0) const;

    //! Prints a keyline to standard out (static)
    /*!
     * Writes
     *         "keyLineTok" = "keyArgTok"
     *
     * @param keyLineTok TOKEN Ptr representing the keyline
     * @param keyArgTok  TOKEN Ptr representing the arguments 
     */
    static void print_keyLine(const TK_TOKEN *keyLineTok, 
			      const TK_TOKEN *keyArgTok);

    //! Indents a line by a certain amount
    /*!
     *  Writes spaces to stdout
     *  @param ilvl number of indent lvls.
     */
    static void print_indent(int ilvl);


    //! This does a recursive search for a Block Entry 
    //! name under the current block
    //! and under all subblocks of the current block.
    /*!
     * It uses a TOKEN ptr as the search input
     *
     * @param nameBN         block name
     * @param includedMatch  Bookean If true then the TOKEN must only be included in the
     *                       Block name. If false, then an exact match is required.
     * @param contribIndex   If the contribIndex is 0 or pos, an exact match with the multiContribIndex() value is
     *                       required. if contribIndex is < 0, then any index is allowed to be matched.
     * @param blockArgName   If this is nonNull then the block argument name is needed to be matched.
     */
    std::set<const BlockEntry*> collectBlockEntries(const TK_TOKEN * const nameBN, bool includedMatch = false, 
						    int contribIndex = -1, const TK_TOKEN * const blockArgName = 0) const;

    //! This does a recursive search for a Block Entry 
    //! name under the current block
    //! and under all subblocks of the current block.
    /*!
     * It uses a TOKEN ptr as the search input
     *
     * @param nameBN         block name
     * @param includedMatch  Bookean If true then the TOKEN must only be included in the
     *                       Block name. If false, then an exact match is required.
     * @param contribIndex   If the contribIndex is 0 or pos, an exact match with the multiContribIndex() value is
     *                       required. if contribIndex is < 0, then any index is allowed to be matched.
     * @param blockArgName   If this is nonNull then the block argument name is needed to be matched.
     */
    std::set<const BlockEntry*> collectBlockEntries(const char *cnameBN, bool includedMatch = false, 
						    int contribIndex = -1, const TK_TOKEN * const blockArgName = 0) const;

    //! Set the int  indicating  whether an error is thrown if an unknown
    //! block or line entry is encountered in the input.
    /*!
     * @param bv value to be set
     *    0 = throw an error if unknown entries are encountered, always
     *    1 = throw an error if unknown entries are encountered unless a block is declared as lenient 
     *    2 = Continue, but print a warning statement.
     *    3 = Continue and don't print anything.
     */
    void set_SkipUnknownEntries(int bv, bool recursive = true);

    //! Get the int  indicating  whether an error is thrown if an unknown
    //! block or line entry is encountered in the input.
    /*!
     *  @return returns the current int
     */
    int get_SkipUnknownEntries() const;

  protected:
   
    //! Pointer to the parent block
    /*!
     * The top block will have this be equal to NULL.
     */
    BlockEntry *ParentBlock;

    //! This is a null terminated list of pointers to LineEntry objects
    //! that define the valid phrases that this block understands
    LineEntry **BlockLineInput;

    /*!
     * This is equal to the number of nonnull items in BlockLineInput
     */
    int numLineInput;

    //! This is a null terminated list of pointers to BlockEntry objects
    //! that define the valid subblocks that this block understands
    BlockEntry **SubBlocks;

    //! Number of subblocks in this list.
    /*!
     * This is equal to the number of nonnull items in SubBlocks
     */
    int numSubBlocks;

    //! Skip unknown entries
    /*!
     *     0  NEVER
     *     1  Sometimes -> (Do the default in s_SkipUnknownEntries
     *     2  Yes, but print
     *     3  Yes, without printing
     *
     *  Default is set to 1 (it defaults to the global default)
     */
    int m_SkipUnknownEntries;

    //! Adjustment to the address originally input to memory assignments
    //! used by line elements in this block. Note, this number is the
    //! cumulative adjustment.
    LONG_PTR m_adjustAddress;

    TK_TOKEN m_ArgTok;

   

  };
}
#endif
