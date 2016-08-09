/**
 * @file BE_BlockEntry.cpp
 *   Definitions for the BlockEntry class which is the base class for all block entry classes
 *   (see \ref blockentryModule and class \link BEInput::BlockEntry BlockEntry\endlink)
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_BlockEntry.h"
#include "LE_LineEntry.h"
#include "mdp_allo.h"
#include <set>

#include <string>

#include <cstdlib>

using namespace std;
using namespace TKInput;
using namespace mdpUtil;
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{

const int BE_ANY_INDEX = -7777779;

//==================================================================================================================================
/*
 *  BlockEntry Constructor
 *
 *
 * Input
 *   FILE *ifp_input = current value of the file pointer
 *                     default is stdin
 *   parentBlock_inp = pointer to the parent of the current block
 *                     default is 0, indicating that this block is the
 *                     top main block that has no parent blocks.
 */
BlockEntry::BlockEntry(const char* blockName, int numTimesRequired, BlockEntry* parentBlock_inp) :
    BaseEntry(blockName, numTimesRequired),
    ParentBlock(parentBlock_inp),
    BlockLineInput(0),
    numLineInput(0),
    SubBlocks(0),
    numSubBlocks(0),
    m_SkipUnknownEntries(1),
    m_adjustAddress(0)
{
    /*
     * Check on sizes
     */
    size_t pz = sizeof(void*);
    size_t lpz = sizeof(LONG_PTR);
    if (pz > lpz) {
        printf("INTERNAL ERROR!\n");
        exit(-1);
    }
}

//==================================================================================================================================
/*
 *
 * BlockEntry(const BlockEntry&) :
 *
 * copy constructor: Deep copy of block structure
 */
BlockEntry::BlockEntry(const BlockEntry& b) :
    BaseEntry(b),
    ParentBlock(b.ParentBlock),
    BlockLineInput(0),
    numLineInput(b.numLineInput),
    SubBlocks(0),
    numSubBlocks(b.numSubBlocks),
    m_SkipUnknownEntries(b.m_SkipUnknownEntries),
    m_adjustAddress(b.m_adjustAddress),
    m_ArgTok(b.m_ArgTok)
{
    BlockLineInput = (LineEntry**) mdp_alloc_ptr_1(numLineInput);
    for (int iLE = 0; iLE < numLineInput; iLE++) {
        BlockLineInput[iLE] =
            b.BlockLineInput[iLE]->duplMyselfAsLineEntry();
    }

    SubBlocks = (BlockEntry**) mdp_alloc_ptr_1(numSubBlocks);
    for (int iSB = 0; iSB < numSubBlocks; iSB++) {
        SubBlocks[iSB] =
            b.SubBlocks[iSB]->duplMyselfAsBlockEntry();
    }
}

//==================================================================================================================================
/*
 *
 *  BlockEntry& operator=(const BlockEntry&)
 */
BlockEntry& BlockEntry::operator=(const BlockEntry& b)
{
    if (&b != this) {
        /*
         * Delete what's there (shared with destructor)
         */
        if (SubBlocks) {
            for (int i = 0; i < numSubBlocks; i++) {
                BlockEntry* sb = SubBlocks[i];
                delete sb;
            }
            mdp_safe_free((void**)&SubBlocks);
        }
        /*
         * Call the destructors for all the lineEntries for the current
         * block
         */
        LineEntry* le;
        if (BlockLineInput) {
            for (int i = 0; i < numLineInput; i++) {
                le = BlockLineInput[i];
                delete le;
            }
            mdp_safe_free((void**)&BlockLineInput);
        }

        BaseEntry::operator=(b);
        ParentBlock = b.ParentBlock;
        numLineInput = b.numLineInput;
        numSubBlocks = b.numSubBlocks;
        m_SkipUnknownEntries = b.m_SkipUnknownEntries;
        m_adjustAddress = b.m_adjustAddress;
        m_ArgTok = b.m_ArgTok;

        BlockLineInput = (LineEntry**) mdp_alloc_ptr_1(numLineInput);
        for (int iLE = 0; iLE < numLineInput; iLE++) {
            BlockLineInput[iLE] =
                b.BlockLineInput[iLE]->duplMyselfAsLineEntry();
        }

        SubBlocks = (BlockEntry**) mdp_alloc_ptr_1(numSubBlocks);
        for (int iSB = 0; iSB < numLineInput; iSB++) {
            SubBlocks[iSB] =
                b.SubBlocks[iSB]->duplMyselfAsBlockEntry();
        }
    }
    return *this;
}

//==================================================================================================================================
/*
 *
 * BlockEntry* duplMyselfAsBlockEntry()
 *
 *     Duplicate the block as a base class
 */
BlockEntry* BlockEntry::duplMyselfAsBlockEntry() const
{
    BlockEntry* newBE = new BlockEntry(*this);
    return newBE;
}

//==================================================================================================================================
/*
 *
 * BaseEntry* duplMyselfAsBaseEntry()
 *
 *     Duplicate the block as a base class, BaseEntry
 */
BaseEntry* BlockEntry::duplMyselfAsBaseEntry() const
{
    BlockEntry* newBE = duplMyselfAsBlockEntry();
    return (BaseEntry*) newBE;
}

//==================================================================================================================================
/*
 *  Block_List Destructor (virtual function):
 *       We may be responsible for releasing memory here.
 *       We must also loop over keylists and subblocks and release
 *       memory.
 */
BlockEntry::~BlockEntry(void)
{
#ifdef DEBUG_DESTRUCTOR
    printf("~BlockEntry called for %s\n", EntryName.orig_str);
#endif
    /*
     * Call the destructors for all subblocks of the current block
     */
    BlockEntry* sb;
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            sb = SubBlocks[i];
            delete sb;
        }
        mdp_safe_free((void**)&SubBlocks);
    }
    /*
     * Call the destructors for all the lineEntries for the current
     * block
     */
    LineEntry* le;
    if (BlockLineInput) {
        for (int i = 0; i < numLineInput; i++) {
            le = BlockLineInput[i];
            delete le;
        }
        mdp_safe_free((void**)&BlockLineInput);
    }
}

//==================================================================================================================================
void BlockEntry::clear()
{
    /*
     * Call the destructors for all subblocks of the current block
     */
    BlockEntry* sb;
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            sb = SubBlocks[i];
            delete sb;
        }
        mdp_safe_free((void**)&SubBlocks);
    }
    numSubBlocks = 0;
    /*
     * Call the destructors for all the lineEntries for the current
     * block
     */
    LineEntry* le;
    if (BlockLineInput) {
        for (int i = 0; i < numLineInput; i++) {
            le = BlockLineInput[i];
            delete le;
        }
        mdp_safe_free((void**)&BlockLineInput);
    }
    numLineInput = 0;
    m_adjustAddress = 0;
}

//==================================================================================================================================
/*
 *
 *  read_block handles the I/O of block commands. It is designed so that
 *  it can be called recursively.
 *
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
 *  followed by an optional name of the block, followed by
 *  an optional "=" sign and an argument list.
 *
 *  Within a block there may be keylines. A keyline is a string
 *  of arbitrary size followed by the "=" key, or a newline,
 *  whichever comes first. After the "=" sign is the argument to
 *  the keyline.
 *
 *
 *    Input
 *  --------------
 *
 *    input_file  - File pointer to read additional keylines
 *                  in the recursive calls to read_block
 *                  Default = stdin
 *
 *    startArgPtr - ptr to the token representing the argument
 *                  to the start block call that initiated the block
 *
 *    parentBlock - pointer to the block class that is the parent
 *                  to the current block. (default no parent indicates
 *                  that this is the top block.
 *
 *    Output
 *   ------------
 *   endArgPtr    - pointer to the TOKEN that represents the
 *                  arguments to the end block line for the current
 *                  block. This is returned to the calling block
 *                  and represents a way to communicate with the
 *                  calling block
 *                   default is the null address.
 *
 */
void
BlockEntry::read_block(FILE* input_file, TK_TOKEN* endArgPtr,
                       const TK_TOKEN* startArgPtr, BlockEntry* parentBlock)
{
    BOOLEAN OK;
    struct TOKEN keyLineTok;  /* Structure to hold the key string itself */
    struct TOKEN keyArgTok;   /* Structure to hold the arguments to the key string  */
    struct TOKEN subBlockEndArg;

    /*
     *   Send an initialization event to the current block
     */
    Initialization(input_file, startArgPtr);

    /*
     *  Figure out whether to skip entries.
     */
    int skipUnknownEntries = m_SkipUnknownEntries;
    if (m_SkipUnknownEntries == 1) {
        skipUnknownEntries = BaseEntry::s_SkipUnknownEntries;
    }
    if (BaseEntry::s_SkipUnknownEntries == 0) {
        skipUnknownEntries = 0;
    }
    /*
     *   Read the next line from the stream that isn't totally made
     *   up of comments. Store everything before the equals sign
     *   in keyLineTok, and everything after the equals sign in keyArgTok
     *    - OK will be false when EOF is read.
     */
    while ((OK = get_next_keyLine(input_file, &keyLineTok, &keyArgTok))) {

#ifdef DEBUG
        print_keyLine(&keyLineTok, &keyArgTok);
#endif
        /*
         *    If the first two keyString tokens are "start block"
         *    handle the nested block call.
         */
        if (keyLineTok.ntokes > 2) {
            if (strmatch(keyLineTok.tok_ptrV[0], "start") &&
                    strmatch(keyLineTok.tok_ptrV[1], "block")) {

                /*
                 *  Find the match to the name of the subblock
                 */
                strip_item_from_token(0, &keyLineTok);
                strip_item_from_token(0, &keyLineTok);
                BlockEntry* subBlockPtr = match_block(&keyLineTok);
                if (!subBlockPtr) {
                    if (skipUnknownEntries >= 2) {
                        if (skipUnknownEntries == 2) {
                            printf("read_block WARNING: Unknown block name:\n");
                            printf("   Enclosing block name = %s, ", EntryName.orig_str);
                            printf("   Unknown block name = %s\n", keyLineTok.orig_str);
                        }
                        skip_block(input_file, &keyLineTok, this);
                    } else {
                        print_usage(2);
                        throw BI_UnknownSubBlock("read_block ("+ string(EntryName.orig_str) + ")",
                                                 "ERROR Unrecognized block: " + string(keyLineTok.orig_str));
                    }
                } else {

                    /*
                     * Do initialization of the subblock at the current block lvl
                     */
                    Initialize_SubBlock(subBlockPtr, &keyArgTok);
                    /*
                     * Call this function for the subblock in a recursive, nested
                     * fashion
                     */
                    subBlockPtr->read_block(input_file, &subBlockEndArg, &keyArgTok, this);
                    /*
                     * Finalize the subblock at the current block lvl
                     */
                    Wrapup_SubBlock(subBlockPtr, &subBlockEndArg);
                }
                /*
                 * Skip to get next line
                 */
                continue;
            } else {
                /*
                 *   Check for the end of the current block -> it had better be
                 *   the block we are currently working on, or it is an error
                 *   - Then, send the block signal to the event
                 *   - and,  break so that we can get out of the function
                 */
                if (strmatch(keyLineTok.tok_ptrV[0], "end") &&
                        strmatch(keyLineTok.tok_ptrV[1], "block")) {
                    strip_item_from_token(0, &keyLineTok);
                    strip_item_from_token(0, &keyLineTok);

                    if (toktokmatch(&keyLineTok, &EntryName)) {
                        Wrapup(input_file, &keyArgTok);
                        copyTokStruct(endArgPtr, &keyArgTok);
                        return;
                    } else {
                        throw BI_EndBlockMismatch("read_block",
                                                  keyLineTok.orig_str,
                                                  this->EntryName.orig_str);
                    }
                }
            }
        }
        /*
         *  Find the match to the name of the keyline
         */
        if (keyLineTok.ntokes > 0) {
            LineEntry* curLE = match_keyLine(&keyLineTok);
            if (!curLE) {
                /*
                 * If we can't find a match, we query the static variable
                 * s_SKipUnknownEntries in order to determine what to do
                 * about it. We either skip it, or throw an error condition, or
                     * write a warning.
                 */
                if (skipUnknownEntries >= 2) {
                    if (skipUnknownEntries == 2)  {
                        printf("read_block WARNING: Unknown line entry:\n");
                        printf("   Enclosing block name = %s, ", EntryName.orig_str);
                        printf("   Unknown line entry name = %s\n", keyLineTok.orig_str);
                    }
                    skip_lineEntry(input_file, &keyLineTok, &keyArgTok);
                } else  {
                    throw BI_UnknownKeyLine("BlockEntry::read_block (" + string(EntryName.orig_str) + ")",
                                            "ERROR Unrecognized keyline: " + string(keyLineTok.orig_str));
                }
            } else {
                /*
                 * Process the line entry
                 */
                process_LineEntry(curLE, &keyArgTok);
            }
        }
    }

    /*
     *     When we are here, either a "block end" signal has occurred, or
     *     an EOF has occurred. Check the variable OK to see which has
     *     occurred.
     *     An EOF in the middle of the main block is not an error
     *      -> just send the end of block signal and leave
     */
    if (!OK) {
        if (parentBlock == NULL) {
            Wrapup(input_file, NULL);
        } else {
            fprintf(stderr,"EOF occurred during the middle of subblock\n");
            exit(-1);
        }
    }
}
//==================================================================================================================================
/************************************************************************
 *
 * read_this_block_only()
 *
 *  This function will skip ahead to read a subblock within a block
 *  structure. Actually, it doesn't check the block structure
 *  until it gets within the actual block.
 */
void BlockEntry::
read_this_block_only(FILE* input_file)
{
    BOOLEAN OK;
    struct TOKEN keyLineTok;  /* Structure to hold the key string itself */
    struct TOKEN keyArgTok;   /* Structure to hold the arguments to the
				 key string  */
    struct TOKEN subBlockEndArg;
    bool found_and_processed = false;


    /*
     *   Read the next line from the stream that isn't totally made
     *   up of comments. Store everything before the equals sign
     *   in keyLineTok, and everything after the equals sign in keyArgTok
     *    - OK will be false when EOF is read.
     */
    while ((OK = get_next_keyLine(input_file, &keyLineTok, &keyArgTok))) {
#ifdef DEBUG
        if (Debug_Flag) {
            print_keyLine(&keyLineTok, &keyArgTok);
        }
#endif
        /*
         *    If the first two keyString tokens are "start block"
         *    handle the nested block call.
         */
        if (keyLineTok.ntokes > 2) {
            if (strmatch(keyLineTok.tok_ptrV[0], "start") &&
                    strmatch(keyLineTok.tok_ptrV[1], "block")) {

                /*
                 *  Find the match to the name of the subblock
                 */
                strip_item_from_token(0, &keyLineTok);
                strip_item_from_token(0, &keyLineTok);
                if (toktokmatch(&EntryName, &keyLineTok)) {
                    /*
                     * Call this function for the subblock in a recursive, nested
                     * fashion
                     *
                     */
                    read_block(input_file, &subBlockEndArg, &keyArgTok, this);
                    found_and_processed = true;
                    break;
                }
            }
        }

        /*
         *     When we are here, either a "block end" signal has occurred, or
         *     an EOF has occurred. Check the variable OK to see which has
         *     occurred.
         */
        if (!OK) {
            if (!found_and_processed) {
                fprintf(stderr,"EOF occurred while searching for specific block\n");
                exit(-1);
            }
        }
    }
}

//==================================================================================================================================

/*
 * skip_block:
 *
 *  skip_block skips a block and all embedded blocks. On entry
 *  the file position should be just after the block header of the
 *  block to be skipped.
 *
 *    Input
 *  --------------
 *
 *    input_file  - File pointer to read additional keylines
 *                  in the recursive calls to read_block
 *                  Default = stdin
 *
 *    BlockName  - ptr to the token representing the name
 *                 of the block to skip.
 *
 *    parentBlock - pointer to the block class that is the parent
 *                  to the current block. (default no parent indicates
 *                  that this is the top block.
 *
 */
void BlockEntry::skip_block(FILE* input_file, TK_TOKEN* blockName,
                            BlockEntry* parentBlock)
{
    BOOLEAN OK;
    int blockLevel = 1;
    struct TOKEN keyLineTok;  /* Structure to hold the key string itself */
    struct TOKEN keyArgTok;   /* Structure to hold the arguments to the
				 key string  */
    struct TOKEN subBlockEndArg;
    while ((OK = get_next_keyLine(input_file, &keyLineTok, &keyArgTok))) {
#ifdef DEBUG
        print_keyLine(&keyLineTok, &keyArgTok);
#endif

        if (keyLineTok.ntokes > 2) {
            if (strmatch(keyLineTok.tok_ptrV[0], "start") &&
                    strmatch(keyLineTok.tok_ptrV[1], "block")) {
                blockLevel++;
            } else {
                /*
                 *   Check for the end of the current block -> it had better be
                 *   the block we were told to skip
                 *   - and,  break so that we can get out of the function
                 */
                if (strmatch(keyLineTok.tok_ptrV[0], "end") &&
                        strmatch(keyLineTok.tok_ptrV[1], "block")) {
                    strip_item_from_token(0, &keyLineTok);
                    strip_item_from_token(0, &keyLineTok);
                    blockLevel--;
                    if (blockLevel == 0) {
                        if (toktokmatch(&keyLineTok, blockName)) {
                            return;
                        } else {
                            throw BI_EndBlockMismatch("skip_block",
                                                      keyLineTok.orig_str,
                                                      this->EntryName.orig_str);
                        }
                    }
                }
            }
        }
    }

    /*
     *     When we are here, either a "block end" signal has occurred, or
     *     an EOF has occurred. Check the variable OK to see which has
     *     occurred.
     *     An EOF in the middle of the main block is not an error
     *      -> just send the end of block signal and leave
     */
    if (!OK) {
        if (parentBlock == NULL) {
            Wrapup(input_file, NULL);
        } else {
            fprintf(stderr,"EOF occurred during the middle of subblock\n");
            fprintf(stderr,"  trying to skip block \"%s\" in block \"%s\"\n",
                    blockName->orig_str, EntryName.orig_str);
            exit(-1);
        }
    }

}
//==================================================================================================================================

// Virtual function called at the start of internally processing the block
/*
 *  This function may be used to start the process of setting up
 *  internal data functions when the block is first
 *  called.
 *  This is also where the current block processes the
 *  arguments specified on the START BLOCK line.
 *  Also, free form input is allowed here as well. This routine may
 *  advance the file pointer, and reassign the starting line of the
 *  block accordingly.
 *  The default behavior is listed below.
 *
 *   -  An Error exit will occur if
 *      a blockARgTok string is actually supplied.
 *
 *   -  Here, we also increment the NumProcessedBlocks counter, which
 *      keeps track of the number of times this particular object has
 *      been called for parent block.
 *
 *   -  A dependency check is made to make sure that the
 *      required dependencies for this block have been satisfied.
 *
 * Derived classes may override the default behavior.
 */
void BlockEntry::Initialization(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* idep = EntryDepList[i];
        if (idep->ResultType() == BIDRT_PT_ERROR) {
            if (!idep->checkDependencySatisfied()) {
                throw BI_InputError("BlockEntry::Initialization",
                                    "Dependency is not satisfied");
            }
        }
        if (idep->ResultType() == BIDRT_ANTITHETICAL_ERROR) {
            if (idep->checkDependencySatisfied()) {
                throw BI_InputError("BlockEntry::Initialization",
                                    "Antithetical dependency is satisfied");
            }
        }
    }
    m_numTimesProcessed++;
    if (blockArgTok->ntokes > 0) {
        m_ArgTok = *blockArgTok;
    }
}

//==================================================================================================================================
/*
 *  Wrapup()  (virtual function)
 *
 *  This function may be used to wrap up the setup of the current
 *  block before returning to the parent block.
 *
 *  Also, free form input is allowed here as well. This routine may
 *  advance the file pointer, and reassign the starting line of the
 *  block accordingly.
 *  The default behavior is to do nothing, but error exit if
 *  a blockARgTok string is actually supplied. Derived classes may
 *  override the default behavior.
 */
void BlockEntry::Wrapup(FILE* ifp_input, const TK_TOKEN* blockArgTok)
{
    if (blockArgTok) {
        if (blockArgTok->ntokes > 0) {
            fprintf(stderr,
                    " Wrapup: Error blockArgTok argument: %s, ignored\n",
                    blockArgTok->orig_str);
            exit(-1);
        }
    }
    /*
     * Check for required line Entries
     */
    if (BlockLineInput) {
        for (int i = 0; i < numLineInput; i++) {
            LineEntry* le = BlockLineInput[i];
            if (! le->checkRequirements(true)) {
                string ss = string(EntryName.orig_str);
                throw BI_MissingRequired(string("BlockEntry::Wrapup for block ") + ss,
                                         le->EntryName.orig_str,
                                         le->m_numTimesRequired,
                                         le->m_numTimesProcessed);
            }
        }
    }
    /*
     * Check for required SubBlocks
     */
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* sb = SubBlocks[i];
            if (! sb->checkRequirements(true)) {
                string ss("SubBlock ");
                ss += sb->EntryName.orig_str;
                throw BI_MissingRequired("BlockEntry::Wrapup", ss,
                                         sb->m_numTimesRequired,
                                         sb->m_numTimesProcessed);
            }
        }
    }
}

//==================================================================================================================================
/*
 * adjustAddress:
 *
 *  Note, the input is a incremental adjustment. The cumulative
 *  adjustment is kept in the class object.
 *
 */
void BlockEntry::adjustAddress(LONG_PTR adjustAAA)
{
    m_adjustAddress += adjustAAA;
    if (BlockLineInput) {
        for (int i = 0; i < numLineInput; i++) {
            LineEntry* le = BlockLineInput[i];
            le->adjustAddress(adjustAAA);
        }
    }
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* be = SubBlocks[i];
            be->adjustAddress(adjustAAA);
        }
    }
}

//==================================================================================================================================
void BlockEntry::set_default(double defValue)
{
    throw BI_InputError("BlockEntry::set_default",
                        "base entry called and not overriden");
}
//==================================================================================================================================
/*
 *
 * ZeroLineCount():
 *
 *  This call zeroes the line count information out. this is needed
 *  when BlockEntry is used in a multiple pass format.
 */
void BlockEntry::ZeroLineCount()
{
    zeroTimesCounter();
    if (BlockLineInput) {
        for (int i = 0; i < numLineInput; i++) {
            LineEntry* le = BlockLineInput[i];
            le->zeroTimesCounter();
        }
    }
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* sb = SubBlocks[i];
            sb->ZeroLineCount();
        }
    }
}

//==================================================================================================================================
/*
 * print_keyLine() (static function):
 *
 *   Prints a line of input to standard out.
 *
 */
void BlockEntry::print_keyLine(const TK_TOKEN* keyLineTok,
                               const TK_TOKEN* keyArgTok)
{
    int i;
    printf("\t\"%s\"", keyLineTok->orig_str);
    if (keyArgTok->ntokes > 0) {
        printf(" =");
        for (i = 0; i < (keyArgTok->ntokes); i++) {
            printf(" \"%s\"", keyArgTok->tok_ptrV[i]);
        }
    }
    printf("\n");
}

//==================================================================================================================================
/*
 * print_indent (BlockEntry static function)
 */
void BlockEntry::print_indent(int ilvl)
{
    for (int i = 0; i < ilvl; i++) {
        printf("    ");
    }
}

/*
 *  print_usage():
 *
 *    This routine will print a summary of the permissible line Entries
 *    and subblocks for this particular block
 */
void BlockEntry::print_usage(int indent_lvl) const
{
    /*
     * The null pointer indicates the top block.
     * print a message for the top block
     */
    if (!ParentBlock) {
        print_indent(indent_lvl);
        printf("*********************************");
        printf(" START OF USAGE INFO FOR TOP BLOCK ");
        printf("*********************************\n");
    } else {
        print_indent(indent_lvl);
        printf("START BLOCK %s", EntryName.orig_str);
        /*
         * Any information about arguments to start block would go here
         */
        string sep = " = ";
        if (m_TimesRequiredType == 1) {
            sep = " >= ";
        } else  if (m_TimesRequiredType == 2 || m_TimesRequiredType == 3) {
            sep = " <= ";
        }
        if (m_numTimesRequired >= 1) {
            if (m_TimesRequiredType == 2) {
                printf("           (OPTIONAL BLOCKS %s %d)", sep.c_str(), m_numTimesRequired);
            } else {
                printf("           (REQUIRED BLOCKS %s %d)", sep.c_str(), m_numTimesRequired);
            }
        }
        if (m_numTimesRequired == 0) {
            printf("           (OPTIONAL BLOCK)");
        }
        printf("\n");
    }
    indent_lvl++;
    if (m_numTimesProcessed || m_multiContribIndex) {
        print_indent(indent_lvl);
        printf("=> timesProcessed = %d ", m_numTimesProcessed);
        printf("multiContribIndex = %d ", m_multiContribIndex);
        printf("\n");
    }
    for (int i = 0; i < numLineInput; i++) {
        LineEntry* le = BlockLineInput[i];
        le->print_usage(indent_lvl);
    }
    for (int i = 0; i < numSubBlocks; i++) {
        BlockEntry* sb = SubBlocks[i];
        sb->print_usage(indent_lvl);
    }
    /*
     * The null pointer indicates the top block.
     * print a message for the top block
     */
    indent_lvl--;
    if (!ParentBlock) {
        print_indent(indent_lvl);
        printf("*********************************");
        printf("  END OF USAGE INFO FOR TOP BLOCK  ");
        printf("*********************************\n");
    } else {
        print_indent(indent_lvl);
        printf("END BLOCK %s", EntryName.orig_str);
        /*
         * Any information about arguments  to end block would go here
         */
        printf("\n");
    }
}

//==================================================================================================================================
/*
 *  match_block():
 *
 *    Match the block name as defined by a sequence of tokens.
 *    Return the address of the block entry.
 *
 *    This function will lop off unwanted trailing "block" "end" and "start"
 *    words from the token list.
 *
 *    If the contribIndex is 0 or pos, an exact match with the index is
 *    required. if contribIndex is < 0, the following procedure is used.
 *
 *    The value of the lowest  multiContribIndex value for a block that
 *    hasn't been processed gets selected.
 *
 */
BlockEntry* BlockEntry::match_block(const TK_TOKEN* keyLinePtr, int contribIndex, bool includedMatch) const
{
    BlockEntry* biBest = 0;
    int multiContribIndexBest = -1;
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* bi = SubBlocks[i];
            bool mmm = false;
            if (includedMatch) {
                mmm = toktokincluded(&(bi->EntryName), keyLinePtr);
            } else {
                mmm = toktokmatch(&(bi->EntryName), keyLinePtr);
            }
            if (mmm) {
                int multiContribIndex = bi->multiContribIndex();
                if (contribIndex == multiContribIndex) {
                    return bi;
                }
                if (contribIndex < 0) {
                    int numTimesProc = bi->get_NumTimesProcessed();
                    if (numTimesProc <= 0) {
                        if (multiContribIndexBest == -1) {
                            biBest = bi;
                            multiContribIndexBest =  multiContribIndex;
                        } else {
                            int numTimesProcBest = biBest->get_NumTimesProcessed();
                            if (numTimesProcBest > 0) {
                                biBest = bi;
                                multiContribIndexBest =  multiContribIndex;
                            } else if (multiContribIndex < multiContribIndexBest) {
                                biBest = bi;
                                multiContribIndexBest =  multiContribIndex;


                            } else if (multiContribIndex == multiContribIndexBest) {
                                throw BI_InputError("BlockEntry::match_block",
                                                    "Two blocks are the same:");
                            }
                        }
                    } else {
                        if (!biBest) {
                            biBest = bi;
                            multiContribIndexBest =  multiContribIndex;
                        }
                    }
                }
            }
        }
    }
    return (biBest);
}

//==================================================================================================================================
BlockEntry* BlockEntry::match_block_argName(const TK_TOKEN* keyLinePtr, bool includedMatch,
                                            int contribIndex, const TK_TOKEN* keyArgName) const
{
    if (!keyArgName) {
        return match_block(keyLinePtr, contribIndex, includedMatch);
    }
    bool argMatch = false;
    BlockEntry* biBest = 0;
    int multiContribIndexBest = -1;
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* bi = SubBlocks[i];
            bool mmm = false;
            if (includedMatch) {
                mmm = toktokincluded(&(bi->EntryName), keyLinePtr);
            } else {
                mmm = toktokmatch(&(bi->EntryName), keyLinePtr);
            }
            if (mmm) {
                if (includedMatch) {
                    argMatch = toktokincluded(&(bi->m_ArgTok), keyArgName);
                } else {
                    argMatch = toktokmatch(&(bi->m_ArgTok), keyArgName);
                }
                if (argMatch) {
                    return bi;
                }
                int multiContribIndex = bi->multiContribIndex();
                if (contribIndex >= 0) {
                    if (contribIndex == multiContribIndex) {
                        biBest = bi;
                    }
                }
                if (contribIndex < 0) {
                    int numTimesProc = bi->get_NumTimesProcessed();
                    if (numTimesProc <= 0) {
                        if ((multiContribIndexBest == -1) ||
                                (multiContribIndex < multiContribIndexBest)) {
                            biBest = bi;
                            multiContribIndexBest =  multiContribIndex;
                        } else   if (multiContribIndex == multiContribIndexBest) {
                            throw BI_InputError("BlockEntry::match_block",
                                                "Two blocks are the same:");
                        }
                    } else {
                        if (!biBest) {
                            biBest = bi;
                        }
                    }
                }
            }
        }
    }
    return (biBest);
}
//=============================================================================================================
/*
 *  match_block():
 *
 *    Match the block name as defined by a character string, which
 *    is translated into a token object.
 *    Return the address of the block entry.
 *
 *    This function will lop off unwanted trailing "block" "end" and "start"
 *    words from the token list.
 */
BlockEntry* BlockEntry::match_block(const char* blockstring, int contribIndex, bool includedMatch) const
{
    TOKEN* t1_ptr = new TOKEN(blockstring);
    BlockEntry* ber = match_block(t1_ptr, contribIndex, includedMatch);
    delete(t1_ptr);
    return ber;
}
//=================================================================================================================
/*
 *  Initialize_SubBlock() (virtual function):
 *
 * Do initialization of the subblock at the current block lvl.
 *
 *  Input
 * ------
 *    subBlockPtr = pointer to the subblock object
 *    keyArgTok   = arguments to the subblock call.
 */
void BlockEntry::Initialize_SubBlock(BlockEntry* subBlockPtr, const TK_TOKEN* keyArgTok)
{
    if (keyArgTok) {
        m_ArgTok = *(keyArgTok);
    } else {
        m_ArgTok = TK_TOKEN("");
    }
}
//===================================================================================================================
/*
 *  Wrapup_SubBlock() (virtual function):
 *
 * Do wrapup of the subblock at the current block lvl.
 *
 *  Input
 * ------
 *    subBlockPtr = pointer to the subblock object
 *    subBlockEndPtr  = arguments to the subblock call.
 */
void BlockEntry::Wrapup_SubBlock(BlockEntry* subBlockPtr,
                                 const TK_TOKEN* subBlockEndPtr)
{
}
//==========================================================================================================
/*
 *
 * match_keyLine():
 *
 *    See if a match exists between the character string and the keyLine
 *    structure.
 */
LineEntry* BlockEntry::match_keyLine(const TK_TOKEN* lineEntryTok, int contribIndex) const
{
    LineEntry* leBest = 0;
    int multiContribIndexBest = -1;
    if (BlockLineInput) {
        for (int i = 0; i < numLineInput; i++) {
            LineEntry* le = BlockLineInput[i];
            if (toktokmatch(&(le->EntryName), lineEntryTok)) {
                int multiContribIndex = le->multiContribIndex();
                if (contribIndex == multiContribIndex) {
                    return le;
                }
                if (multiContribIndex > multiContribIndexBest) {
                    leBest = le;
                    multiContribIndexBest = multiContribIndex;
                } else if (multiContribIndex == multiContribIndexBest) {
                    throw BI_InputError("BlockEntry::match_keyLine",
                                        "Two LineEntries are the same:");
                }
            }
        }
    }
    return (leBest);
}
//==========================================================================================================
/*
 * process_LineEntry() (virtual function):
 *
 *   The basic approach is to pass the control of dealing with line entries
 *   down to the individual LineEntry objects themselves.
 */
void BlockEntry::process_LineEntry(LineEntry* curLE,
                                   const TK_TOKEN* keyArgTokPtr)
{
    curLE->process_LineEntry(keyArgTokPtr);
    if (BaseEntry::s_PrintProcessedLine) {
        curLE->print_ProcessedLine(keyArgTokPtr);
    }
}

//==================================================================================================================================
/*
 * skip_lineEntry() (virtual function):
 *
 *  Currently, this is just a holder for perhaps printing more
 *  debugging information.
 */
void BlockEntry::skip_lineEntry(FILE* ifp, const TK_TOKEN* keyLineTok,
                                const TK_TOKEN* keyArgTokPtr) const
{
}

//==================================================================================================================================
/*
 * addSubBlock
 *
 *   This adds a subblock entry structure to the current BlockEntry
 *   structure.
 */
void BlockEntry::addSubBlock(BlockEntry* sb)
{
    if (!sb) {
        return;
    }
    /*
     * Search for an existing match. If found delete the
     * old entry and replace with new
     */
    int multiContribIndexNew = sb->multiContribIndex();
    for (int iLE = 0; iLE < numSubBlocks; iLE++) {
        BlockEntry* sbOld = SubBlocks[iLE];
        if (toktokmatch(&(sb->EntryName), &(sbOld->EntryName))) {
            int multiContribIndexOld = sbOld->multiContribIndex();
            if (multiContribIndexNew == multiContribIndexOld) {
#ifdef DEBUG_HKM
                TOKEN* sbt = &(sb->EntryName);
                printf("ERROR addSubBlock: overwriting an existing subblock: "
                       "\"%s\"\n", sbt->orig_str);
                exit(-1);
#endif
                delete sbOld;
                SubBlocks[iLE] = sb;
                return;
            }
        }
    }
    mdp_realloc_ptr_1((void***)&SubBlocks, numSubBlocks+2, numSubBlocks);
    SubBlocks[numSubBlocks] = sb;
    SubBlocks[numSubBlocks+1] = 0;
    numSubBlocks++;
    sb->ParentBlock = this;
}

//==================================================================================================================================
/*
 * addLineEntry
 *
 *   This adds a Line Entry structure to the current BlockEntry
 *   structure.
 */
void BlockEntry::addLineEntry(LineEntry* le)
{
    if (!le) {
        return;
    }
    /*
     * Search for an existing match. If found delete the
     * old entry and replace with new
     */
    int multiContribIndexNew = le->multiContribIndex();
    for (int iLE = 0; iLE < numLineInput; iLE++) {
        LineEntry* leOld = BlockLineInput[iLE];
        if (toktokmatch(&(le->EntryName), &(leOld->EntryName))) {
            int multiContribIndexOld = leOld->multiContribIndex();
            if (multiContribIndexNew == multiContribIndexOld) {
                if (le == leOld) {
                    throw BI_InputError("BlockEntry::addLineEntry",
                                        "Adding same line entry to object again:" +
                                        string(le->EntryName.orig_str));
                }
                delete leOld;
                BlockLineInput[iLE] = le;
                return;
            }
        }
    }
    /*
     * No match was found. Return
     */
    mdp_realloc_ptr_1((void***)&BlockLineInput,
                      numLineInput + 2, numLineInput);
    BlockLineInput[numLineInput] = le;
    BlockLineInput[numLineInput+1] = 0;
    numLineInput++;
}

//==================================================================================================================================
/*
 *  reportNumProcessedLines
 *
 * This does a recursive search for a Line Entry on the current block
 * and all subblocks of the current block, using a character string
 * as the input.
 */
int BlockEntry::reportNumProcessedLines(const char* lineName) const
{
    if (!lineName) {
        return 0;
    }
    LineEntry* le = searchLineEntry(lineName);
    if (le) {
        return le->m_numTimesProcessed;
    }
    return 0;
}

//==================================================================================================================================
/*
 *
 * checkRequirements:
 *
 * Check for all requirements being met at the end of input.
 *
 * @param throwSpecificError If true then you should throw a
 *        specific error condition, if you have one. If not,
 *        an generic error condition will be thrown on return.
 */
bool BlockEntry::checkRequirements(bool throwSpecificError)
{
    /*
     * The first thing we do is to check for satisfaction
     * of any Dependency Result Types that might zero
     * out the Number of times required field.
     * Thus, this may have been a required block. However,
     * after processing another block or line, this block
     * is deemed to no longer be needed, or to be an
     * optional block.
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* dep = EntryDepList[i];
        if (dep->ResultType() == BIDRT_ZERONUMTIMESREQUIRED) {
            if (dep->checkDependencySatisfied()) {
                m_numTimesRequired = 0;
            }
        }
    }
    /*
     * The next thing we do is to check for satifaction
     * of any dependency result types that might turn on
     * the number of times required field for this block.
     * Thus, if this may have been an optional block. However,
     * after a card has been introduced elsewhere, this
     * block is now a required block.
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* dep = EntryDepList[i];
        if (dep->ResultType() == BIDRT_ONENUMTR) {
            if (dep->checkDependencySatisfied()) {
                m_numTimesRequired = 1;
            }
        }
    }

    /*
     * The next thing we do is to check for satisfaction
     * of any Dependency Result Types that might
     * be antithetical to the current card.
     */
    for (int i = 0; i < NumEntryDependencies; i++) {
        BI_Dependency* dep = EntryDepList[i];
        if (dep->ResultType() == BIDRT_ANTITHETICAL_ERROR) {
            if (dep->checkDependencySatisfied()) {
                throw BI_InputError("BlockEntry::checkRequirements()",
                                    "Mutually exclusive BaseEntry objects were invoked");
            }
        }
    }

    if (m_numTimesRequired) {
        if (m_TimesRequiredType == 0) {
            if (m_numTimesRequired != m_numTimesProcessed) {
                return false;
            }
        } else if (m_TimesRequiredType == 1) {
            if (m_numTimesRequired >  m_numTimesProcessed) {
                return false;
            }
        } else if (m_TimesRequiredType == 2) {
            if (m_numTimesRequired < m_numTimesProcessed) {
                return false;
            }
        } else if (m_TimesRequiredType == 3) {
            if (m_numTimesRequired < m_numTimesProcessed) {
                return false;
            }
            if (m_numTimesProcessed == 0) {
                return false;
            }
        }
    }
    return true;
}

//==================================================================================================================================
/*
 * searchLineEntry:
 *
 *    This does a recursive search for a LineEntry on the current block
 * and all subblocks of the current block, using a character string
 * as the input.
 */
LineEntry* BlockEntry::
searchLineEntry(const char* const lineName) const
{
    if (!lineName) {
        return 0;
    }
    TOKEN* tok_ptr = new TOKEN(lineName);
    LineEntry* le = searchLineEntry(tok_ptr);
    delete tok_ptr;
    return le;
}

//==================================================================================================================================
/*
 * searchLineEntry:
 *
 *  This does a recursive search for a LineEntry on the current block
 * and all subblocks of the current block.
 */
LineEntry* BlockEntry::searchLineEntry(const TK_TOKEN* const nameLE) const
{
    LineEntry* le = NULL;
    if (BlockLineInput) {
        le = match_keyLine(nameLE);
        if (le) {
            return (le);
        }
    }
    if (SubBlocks) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* bi = SubBlocks[i];
            le = bi->searchLineEntry(nameLE);
            if (le) {
                return (le);
            }
        }
    }
    return (NULL);
}
//==================================================================================================================
/*
 * searchBlockEntry:
 *
 * This does a recursive search for a Block Entry
 * name under the current block
 * and under all subblocks of the current block.
 */
BlockEntry* BlockEntry::searchBlockEntry(const char* const bName, bool includedMatch,
                                         int contribIndex, const TK_TOKEN* const blockArgName) const
{
    if (!bName) {
        return 0;
    }
    TOKEN* tok_ptr = new TOKEN(bName);
    BlockEntry* be = searchBlockEntry(tok_ptr, includedMatch, contribIndex, blockArgName);
    delete tok_ptr;
    return be;
}
//==================================================================================================================================
/*
 * searchBlockEntry:
 *
 * This does a recursive search for a Block Entry
 * name under the current block
 * and under all subblocks of the current block.
 */
BlockEntry* BlockEntry::searchBlockEntry(const TK_TOKEN* const nameBN, bool includedMatch,
                                         int contribIndex, const TK_TOKEN* const blockArgName) const
{
    BlockEntry* bi = NULL;
    if (includedMatch) {
        if (toktokincluded(&EntryName, nameBN)) {
            if (blockArgName) {
                if (toktokincluded(&m_ArgTok, blockArgName)) {
                    return (const_cast<BlockEntry*>(this));
                }
            } else {
                if (contribIndex >= 0) {
                    if (contribIndex == bi->multiContribIndex()) {
                        return (const_cast<BlockEntry*>(this));
                    }
                } else {
                    return (const_cast<BlockEntry*>(this));
                }
            }
        }
        if (SubBlocks) {
            bi = match_block_argName(nameBN,  true, contribIndex, blockArgName);
            if (bi) {
                return bi;
            }
            for (int i = 0; i < numSubBlocks; i++) {
                BlockEntry* sbi = SubBlocks[i];
                bi = sbi->searchBlockEntry(nameBN, true, contribIndex, blockArgName);
                if (bi) {
                    return (bi);
                }
            }
        }
    } else {
        if (toktokmatch(&EntryName, nameBN)) {
            if (blockArgName) {
                if (toktokincluded(&m_ArgTok, blockArgName)) {
                    return (const_cast<BlockEntry*>(this));
                }
            } else {
                if (contribIndex >= 0) {
                    if (contribIndex == bi->multiContribIndex()) {
                        return (const_cast<BlockEntry*>(this));
                    }
                } else {
                    return (const_cast<BlockEntry*>(this));
                }
            }
        }

        if (SubBlocks) {
            bi = match_block_argName(nameBN, false, contribIndex, blockArgName);
            if (bi) {
                return bi;
            }
            for (int i = 0; i < numSubBlocks; i++) {
                BlockEntry* sbi = SubBlocks[i];
                bi = sbi->searchBlockEntry(nameBN, false, contribIndex, blockArgName);
                if (bi) {
                    return (bi);
                }
            }
        }
    }
    return (NULL);
}
//=================================================================================================================
std::set<const BlockEntry*> BlockEntry::collectBlockEntries(const char* cnameBN, bool includedMatch,
                                                            int contribIndex, const TK_TOKEN* const blockArgName) const
{
    std::set<const BlockEntry*> cc;
    if (!cnameBN) {
        return cc;
    }
    TOKEN tok_ptr(cnameBN);
    cc = collectBlockEntries(&tok_ptr, includedMatch, contribIndex,  blockArgName);
    return cc;
}
//=================================================================================================================

// This does a recursive search for a Block Entry name
// and under all subblocks of the current block.
/*
 * It uses a TOKEN ptr as the search input
 *
 * @param nameBN         block name
 * @param includedMatch  Bookean If true then the TOKEN must only be included in the
 *                       Block name. If false, then an exact match is required.
 * @param contribIndex   If the contribIndex is 0 or pos, an exact match with the multiContribIndex() value is
 *                       required. if contribIndex is < 0, then any index is allowed to be matched.
 * @param blockArgName   If this is nonNull then the block argument name is needed to be matched.
 */
std::set<const BlockEntry*> BlockEntry::collectBlockEntries(const TK_TOKEN* const nameBN, bool includedMatch,
                                                            int contribIndex, const TK_TOKEN* const blockArgName) const
{
    std::set<const BlockEntry*> cc;
    std::set<const BlockEntry*> cc_sub;

    if (includedMatch) {
        if (toktokincluded(&EntryName, nameBN)) {
            if (blockArgName) {
                if (toktokincluded(&m_ArgTok, blockArgName)) {
                    cc.insert(const_cast<BlockEntry*>(this));
                }
            } else {
                if (contribIndex >= 0) {
                    if (contribIndex == multiContribIndex()) {
                        cc.insert(const_cast<BlockEntry*>(this));
                    }
                } else {
                    cc.insert(const_cast<BlockEntry*>(this));
                }
            }
        }
        if (SubBlocks) {
            for (int i = 0; i < numSubBlocks; i++) {
                BlockEntry* sbi = SubBlocks[i];
                cc_sub = sbi->collectBlockEntries(nameBN, true, contribIndex, blockArgName);
                cc.insert(cc_sub.begin(), cc_sub.end());
            }
        }
    } else {
        if (toktokmatch(&EntryName, nameBN)) {
            if (blockArgName) {
                if (toktokincluded(&m_ArgTok, blockArgName)) {
                    cc.insert(const_cast<BlockEntry*>(this));
                }
            } else {
                if (contribIndex >= 0) {
                    if (contribIndex == multiContribIndex()) {
                        cc.insert(const_cast<BlockEntry*>(this));
                    }
                } else {
                    cc.insert(const_cast<BlockEntry*>(this));
                }
            }
        }

        if (SubBlocks) {
            for (int i = 0; i < numSubBlocks; i++) {
                BlockEntry* sbi = SubBlocks[i];
                cc_sub = sbi->collectBlockEntries(nameBN, true, contribIndex, blockArgName);
                if (cc_sub.size() > 0) {
                    cc.insert(cc_sub.begin(), cc_sub.end());
                }
            }
        }
    }
    return (cc);
}
//=================================================================================================================
void BlockEntry::set_SkipUnknownEntries(int bv, bool recursive)
{
    m_SkipUnknownEntries = bv;
    if (recursive) {
        for (int i = 0; i < numSubBlocks; i++) {
            BlockEntry* sbi = SubBlocks[i];
            sbi->set_SkipUnknownEntries(bv, recursive);
        }
    }
}
//=================================================================================================================
int  BlockEntry::get_SkipUnknownEntries() const
{
    return m_SkipUnknownEntries;
}
//=================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
