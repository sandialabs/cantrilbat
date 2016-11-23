/**
 * @file BI_InputError.cpp
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#include "BI_InputError.h"
#include "BaseEntry.h"


#include <cstdlib>
#include <cstring>

using namespace TKInput;

namespace BEInput
{
//==================================================================================================================================
//!  Given an integer, this function returns a C++ string for the integer
/*!
 *       @param[in]   n      Integer to be converted
 *
 *       @return             Converted integer as a string
 */
static std::string integerToStr(int n)
{
    char buf[30];
    (void) sprintf(buf, "%d", n);
    return std::string(buf);
}
//==================================================================================================================================
BI_InputError::BI_InputError(std::string procedure, std::string msg) :
    m_procedure_msg(procedure + ": " + msg)
{
}
//==================================================================================================================================
BI_InputError::~BI_InputError() throw() 
{
}
//==================================================================================================================================
std::string BI_InputError::errorMessage() const
{
    return m_procedure_msg;
}
//==================================================================================================================================
void BI_InputError::append(std::string msg)
{
    m_procedure_msg += msg;
}
//==================================================================================================================================
const char* BI_InputError::what() const throw()
{
    return m_procedure_msg.c_str();
}
//==================================================================================================================================
BI_UnknownKeyLine::BI_UnknownKeyLine(std::string procedure, std::string lineTok) :
    BI_InputError(procedure,
                  "\tUnknown lineKey ->\"" + lineTok + "\"<-")
{
}
//==================================================================================================================================
BI_UnknownSubBlock::BI_UnknownSubBlock(std::string procedure,
                                       std::string lineTok) :
    BI_InputError(procedure,
                  "\tUnknown SubBlock ->\"" + lineTok + "\"<-")
{
}
//==================================================================================================================================
BI_UnknownListEntry::BI_UnknownListEntry(std::string procedure,
                                         std::string lineTok) :
    BI_InputError(procedure, "\tUnknown ListEntry ->\"" + lineTok + "\"<-")
{
}
//==================================================================================================================================
BI_EndBlockMismatch::BI_EndBlockMismatch(std::string procedure,
                                         std::string lineTok,
                                         std::string expected) :
    BI_InputError(procedure, "\tEnd Block Mismatch: got ->\"" + lineTok + "\"<-"
                  " Expected -> \"" + expected + "\"<-")
{
}
//==================================================================================================================================
BI_MissingRequired::BI_MissingRequired(std::string procedure,
                                       std::string lineTok,
                                       int numR, int numP) :
    BI_InputError(procedure, "\tMissing Required Number of Line Entries ->\""
                  + lineTok + "\"<- : Required = " + integerToStr(numR) + " : Processed = " + integerToStr(numP))
{
}
//==================================================================================================================================
BI_MissingRequiredVec::
BI_MissingRequiredVec(std::string procedure, std::string lineTok, int numVneeded, int numVfound) :
    BI_InputError(procedure, "\tMissing Required Number of vector entries ->\""
                  + lineTok + "\"<- :" + integerToStr(numVneeded) + " : " + integerToStr(numVfound))
{
}
//==================================================================================================================================
void BI_SetPrintLevel(int printLevel)
{
    if (printLevel < 2) {
        BaseEntry::set_printProcessedLine(false);
    } else {
        BaseEntry::set_printProcessedLine(true);
    }
    if (printLevel < 1) {
        set_tok_input_print_flag(0);
    } else {
        set_tok_input_print_flag(1);
    }
}
//==================================================================================================================================
}
