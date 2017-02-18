/**
 * @file m1d_globals.cpp
 * This file contains globally defined functions and data
 * that are not part of any class.
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_globals.h"
#include "m1d_app.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdio>
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

int s_printLvl_TimingInformation = 1;

int s_printLvl_DebugTables = 1;
//======================================================================================================================
// By default predictor_corrector printing is turned on, at least to the printLvl_ level.
//==================================================================================================================================
void readEnvironmentalVariables()
{
    char* PC_PRINTING = getenv("M1D_PRINT_TIMINGINFORMATION");
    if (PC_PRINTING) {
        if (PC_PRINTING[0] == 'f' || PC_PRINTING[0] == 'F' || PC_PRINTING[0] == '0') {
            m1d::s_printLvl_TimingInformation = 0;
        } else  if (PC_PRINTING[0] == 't' || PC_PRINTING[0] == 'T' || PC_PRINTING[0] == '1') {
            m1d::s_printLvl_TimingInformation = 1;
        }
    }
    PC_PRINTING = getenv("M1D_PRINT_DEBUGTABLES");
    if (PC_PRINTING) {
        if (PC_PRINTING[0] == 'f' || PC_PRINTING[0] == 'F' || PC_PRINTING[0] == '0') {
            m1d::s_printLvl_DebugTables = 0;
        } else  if (PC_PRINTING[0] == 't' || PC_PRINTING[0] == 'T' || PC_PRINTING[0] == '1') {
            m1d::s_printLvl_DebugTables = 1;
        }
    }
}
//==================================================================================================================================
Appl* app()
{
    return Appl::Instance();
}
//==================================================================================================================================
void
setError(std::string r, std::string msg)
{
    app()->addError(r, msg);
}
//==================================================================================================================================
void
showErrors(std::ostream& f)
{
    app()->showErrors(f);
}
//==================================================================================================================================
void
popError()
{
    app()->popError();
}
//==================================================================================================================================
std::string
intToString(const int p)
{
    char buf[64];
    sprintf(buf, "%d", p);
    return std::string(buf);
}
//==================================================================================================================================
bool checkDblAgree(double d1, double d2, double rtol, double atol)
{
    double dfabs = fabs(d1 - d2);
    double denom = 0.5 * (fabs(d1) + fabs(d2)) * rtol;
    denom = std::max(denom, atol);
    if (dfabs < denom) {
        return true;
    }
    return false;
}
//==================================================================================================================================
void
fwriteTecplotVector(FILE* ofp, const std::vector<double>& vals, int numD, int numPerLine)
{
    //static const char* form = "%20.13E ";
    //const int numD = 13;
    char form[16];
    sprintf(form, "%%%u.%uE ", 7+numD, numD);
    int rsize = 0;
    for (size_t i = 0; i < vals.size(); ++i) {
        fprintf(ofp, form, vals[i]);
        if (++rsize >= numPerLine) {
            fprintf(ofp, "\n");
            rsize = 0;
        }
    }
    if (rsize != 0) {
        fprintf(ofp, "\n");
    }
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

