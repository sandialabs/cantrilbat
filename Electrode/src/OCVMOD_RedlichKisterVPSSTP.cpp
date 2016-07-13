/**
 *  @file RedlichKisterVPSSTP.cpp
 *   Definitions for ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations related to RedlichKister
 *   expansions (see \ref thermoprops
 *    and class \link Zuzax::RedlichKisterVPSSTP RedlichKisterVPSSTP\endlink).
 *
 */
/*
 * Copyright (2009) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "OCVMOD_RedlichKisterVPSSTP.h"


using namespace std;

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//========================================================================================================================
OCVMOD_RedlichKisterVPSSTP::OCVMOD_RedlichKisterVPSSTP() :
    RedlichKisterVPSSTP(),
    rsd_(0)
{
}

//========================================================================================================================
OCVMOD_RedlichKisterVPSSTP::OCVMOD_RedlichKisterVPSSTP(const std::string& inputFile, const std::string& id) :
    RedlichKisterVPSSTP(inputFile, id),
    rsd_(0)
{
}

//========================================================================================================================
OCVMOD_RedlichKisterVPSSTP::OCVMOD_RedlichKisterVPSSTP(XML_Node& phaseRoot, const std::string& id) :
    RedlichKisterVPSSTP(phaseRoot, id),
    rsd_(0)
{
}

//========================================================================================================================
OCVMOD_RedlichKisterVPSSTP::OCVMOD_RedlichKisterVPSSTP(const OCVMOD_RedlichKisterVPSSTP& b) :
    RedlichKisterVPSSTP(),
    rsd_(0)
{
    operator=(b);
}
//========================================================================================================================
OCVMOD_RedlichKisterVPSSTP& OCVMOD_RedlichKisterVPSSTP::operator=(const OCVMOD_RedlichKisterVPSSTP& b)
{
    if (&b == this) {
        return *this;
    }
    OCVMOD_RedlichKisterVPSSTP::operator=(b);
    rsd_ = b.rsd_;
    return *this;
}
//========================================================================================================================
void OCVMOD_RedlichKisterVPSSTP::getChemPotentials(double * mu) const
{
    deriveGibbsCorrection();
    std::copy(m_mu.begin(), m_mu.end(), mu);
}
//=======================================================================================================================
void OCVMOD_RedlichKisterVPSSTP::setup_Start(std::vector<double>& stoichVector, RSD_OCVmodel* OCVmodel)
{
    OCVmodel_ = OCVmodel;
     stoichRxn_ =  stoichVector;
}
//========================================================================================================================
// usually this is Lithium metal only
void OCVMOD_RedlichKisterVPSSTP::setup_AddThermoPhase(std::vector<double>& stoichVector, ThermoPhase *tp)
{
   extraTPList_.push_back(tp);
   size_t nsp = tp->nSpecies();
   std::vector<double> mu(nsp, 0.0);
   m_muExtraList_.push_back(mu);
   extraStoichRxn_.push_back(stoichVector);
}
//========================================================================================================================
void OCVMOD_RedlichKisterVPSSTP::deriveGibbsCorrection() const
{
    ThermoPhase* tp;
    size_t numExtraPhases = extraTPList_.size();
    for (size_t i = 0; i < numExtraPhases; i++) {
         tp = extraTPList_[i];
         std::vector<double>& m_muExtra_phase = m_muExtraList_[i];
         tp->getChemPotentials(DATA_PTR(m_muExtra_phase)); 
    }
}
//========================================================================================================================

} // end of namespace
