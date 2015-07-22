/**
 * @file BE_UnitConversion.cpp
 *   Definitions for the base object that handles unit conversions 
 *  ((see \ref blockentryModule and class \link BEInput::BE_UnitConversion BE_UnitConversion\endlink).
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "BE_UnitConversion.h"
#include "BI_InputError.h"
#include <new>
using std::string;
//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
//==================================================================================================================================
//! Avogadro's number in molecules per kmol
static const double Avogadro =  6.022136736e26;
//==================================================================================================================================
/*
 * BE_UnitConversion() Constructor:
 * This general map will convert to SI units
 * in most cases.
 */
BE_UnitConversion::BE_UnitConversion()
{
    // length
    m_u["m"]        = 1.0;
    m_u["cm"]       = 0.01;
    m_u["km"]       = 1.0e3;
    m_u["mm"]       = 1.0e-3;
    m_u["micron"]   = 1.0e-6;
    m_u["nm"]       = 1.0e-9;
    m_u["A"]        = 1.0e-10;
    m_u["inches"]   = 0.0254;
    m_u["ft"]       = 0.3048;

    // energy
    m_u["J"]        = 1.0;
    m_u["kJ"]       = 1.0e3;
    m_u["cal"]      = 4.184;
    m_u["kcal"]     = 4184.0;
    m_u["eV"]       = 1.602e-19;

    // electricity
    m_u["Amp"]        = 1.0;
    m_u["mAmp"]       = 1.0e-3;
    m_u["V"]        = 1.0;
    m_u["mV"]       = 1.0e-3;

    // quantity
    m_u["mol"]      = 1.0e-3;
    m_u["gmol"]     = 1.0e-3;
    m_u["mole"]     = 1.0e-3;
    m_u["kmol"]     = 1.0;
    m_u["molec"]    = 1.0/Avogadro;

    // temperature -> converter only does multiplication
    m_u["K"]        = 1.0;
    m_u["C"]        = 1.0;

    // mass
    m_u["g"]        = 1.0e-3;
    m_u["kg"]       = 1.0;

    // pressure
    m_u["atm"]      = 1.01325e5;
    m_u["bar"]      = 1.0e5;
    m_u["Pa"]       = 1.0;
    m_u["Pascal"]   = 1.0;
    m_u["torr"]     = 1.01325E5 / 760.;

    // time
    m_u["s"]        = 1.0;
    m_u["min"]      = 60.0;
    m_u["hr"]       = 3600.0;
    m_u["ms"]       = 0.001;
}
//==================================================================================================================================
BE_UnitConversion::BE_UnitConversion(const BE_UnitConversion& right) :
    m_u(right.m_u)
{
}
//==================================================================================================================================
BE_UnitConversion& BE_UnitConversion::operator=(const BE_UnitConversion& right)
{
    if (this != &right) {
        m_u = right.m_u;
    }
    return *this;
}
//==================================================================================================================================
BE_UnitConversion* BE_UnitConversion::duplMyselfAsUnitConversion() const
{
    BE_UnitConversion* bec = new BE_UnitConversion(*this);
    return bec;
}
//==================================================================================================================================
/*
 * The function converts a string expression to
 * a conversion value. The following syntax is used
 * / -> stands for division
 * - -> stands for multiplication
 * [1-9] -> stands for the powers of the previous string
 *          expression.
 */
double BE_UnitConversion::toSI(std::string unitString) const
{
    if (unitString == "") {
        return 1.0;
    }
    double f = 1.0, fctr = 0.0;
    int tsize;
    string u = unitString, tok, tsub;
    string::size_type k;
    char action = '-';
    while (1 > 0) {
        k = u.find_first_of("/-*");
        if (k != string::npos) {
            tok = u.substr(0,k);
        } else {
            tok = u;
        }
        tsize = static_cast<int>(tok.size());
        if (tsize == 0) {
            fctr = 1.0;
        } else if (tok[tsize - 1] == '2') {
            tsub = tok.substr(0,tsize-1);
            fctr = m_u[tsub];
            fctr *= fctr;
        } else if (tok[tsize - 1] == '3') {
            tsub = tok.substr(0,tsize-1);
            fctr = m_u[tsub];
            fctr *= fctr*fctr;
        } else if (tok[tsize - 1] == '4') {
            tsub = tok.substr(0,tsize-1);
            fctr = m_u[tsub];
            fctr *= fctr*fctr*fctr;
        } else if (tok[tsize - 1] == '5') {
            tsub = tok.substr(0,tsize-1);
            fctr = m_u[tsub];
            fctr *= fctr*fctr*fctr*fctr;
        } else if (tok[tsize - 1] == '6') {
            tsub = tok.substr(0,tsize-1);
            fctr = m_u[tsub];
            fctr *= fctr*fctr*fctr*fctr*fctr;
        } else {
            tsub = tok;
            fctr = m_u[tok];
        }

        if (fctr == 0.0) {
            throw BI_InputError("BE_UnitConversion::toSI", "unknown unit: " + tsub);
        }
        if (action == '-' || action == '*') {
            f *= fctr;
        } else if (action == '/') {
            f /= fctr;
        }
        if (k == string::npos) {
            break;
        }
        action = u[k];
        u = u.substr(k+1, u.size());
    }
    return f;
}
//==================================================================================================================================
std::string BE_UnitConversion::returnUsage() const
{
    string hhh;
    hhh = "Units conversion to mks on a per unit_name basis using #/- operators";
    hhh += " (Default: mks units = 1.0)";
    return hhh;
}
//==================================================================================================================================
BE_UnitConversion::~BE_UnitConversion()
{
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
