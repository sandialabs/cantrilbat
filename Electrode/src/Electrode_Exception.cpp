/**
 * @file Electrode_Exception.cpp
 *
 */


#include "Electrode_Exception.h"

namespace Cantera
{

//==============================================================================
Electrode_Error::Electrode_Error(const std::string &proc, const std::string &msg) :
    CanteraError("Electrode_Error: " + proc, msg)
{
}
//==============================================================================
Electrode_Error::~Electrode_Error() throw()
{
}
//==============================================================================
Electrode_Error::Electrode_Error() :
    CanteraError()
{
}
//==============================================================================
}
