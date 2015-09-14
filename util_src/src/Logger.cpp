/**
 * @file Logger.cpp
 * Header for Base class for 'loggers' that write text messages to log files
 */

#include "Logger.h"

namespace mdpUtil {
//----------------------------------------------------------------------------------------------------------------------------------
//==================================================================================================================================
Logger::Logger() :
   out_(&std::cout)
{
}
//==================================================================================================================================
Logger::~Logger()
{
}
//==================================================================================================================================
void Logger::setOutputStream(std::ostream& outputLocation)
{
    out_ = &outputLocation;
}
//==================================================================================================================================
void Logger::write(const std::string& msg)
{
    (*out_) << msg;
}
//==================================================================================================================================
void Logger::writeendl()
{
    (*out_) << std::endl;
}
//==================================================================================================================================
void Logger::error(const std::string& msg)
{
    (*out_) << msg << std::endl;
    exit(-1);
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------