/**
 *  @file cttables_kin.h
 *
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CTTABLES_KIN_H
#define CTTABLES_KIN_H

#include "cttables.h"
#include "zuzax/kinetics.h"
#include "zuzax/kinetics/ExtraGlobalRxn.h"
#include "zuzax/numerics/DenseMatrix.h"

//==================================================================================================================================
class RxnMolChangeLocal
{
public:
    RxnMolChangeLocal(Zuzax::Kinetics* kinPtr, int irxn);
    virtual ~RxnMolChangeLocal();
    RxnMolChangeLocal(Zuzax::Kinetics* kinPtr, Zuzax::ExtraGlobalRxn* egr);

    std::vector<double> m_phaseMoleChange;
    std::vector<double> m_phaseMassChange;
    std::vector<double> m_phaseChargeChange;
    std::vector<double> m_phasePotentials;
    std::vector<int>    m_phaseTypes;
    std::vector<int>    m_phaseDims;

    int m_nPhases;
    Zuzax::Kinetics* m_kinBase;
    int m_iRxn;
    double m_ChargeTransferInRxn;
    double m_beta;
    Zuzax::ExtraGlobalRxn* m_egr;
};

//==================================================================================================================================
class RxnTempTableStuff
{
public:
    RxnTempTableStuff(int irxn, int irxnGE = -1);
    ~RxnTempTableStuff();
    std::vector<double> kfwd_Table;
    std::vector<double> krev_Table;
    std::vector<double> deltaGss_Table;
    std::vector<double> deltaHss_Table;
    std::vector<double> deltaSss_Table;
    std::vector<double> deltaG_Table;
    std::vector<double> deltaH_Table;
    std::vector<double> deltaS_Table;
    std::vector<double> Afwd_Table;
    std::vector<double> EafwddivR_Table;
    std::vector<double> Arev_Table;
    std::vector<double> EarevdivR_Table;
    std::vector<double> kfwdPrime_Table;
    std::vector<double> krevPrime_Table;

    std::vector<double> NetROP_Table;
    std::vector<double> FwdROP_Table;
    std::vector<double> RevROP_Table;
    std::vector<double> Anet_Table;
    std::vector<double> EanetdivR_Table;
    int m_irxn;
    int m_irxnGlobalExtra;
};

//==================================================================================================================================

class TemperatureTable;

//==================================================================================================================================
void doKineticsTablesHomog(Zuzax::PhaseList* pl, Zuzax::Kinetics* gKinetics, TemperatureTable& TT);

void doKineticsTablesHetero(Zuzax::PhaseList* pl, Zuzax::InterfaceKinetics* gKinetics, TemperatureTable& TT);

void processCurrentVsPotTable(Zuzax::RxnMolChange* rmc,
                              Zuzax::PhaseList* pl,
                              int irxn,
                              TemperatureTable& TT,
                              Zuzax::Kinetics& kin,
                              Zuzax::DenseMatrix& kfwd_Table,
                              Zuzax::DenseMatrix& krev_Table,
                              Zuzax::DenseMatrix& deltaG_Table,
                              Zuzax::DenseMatrix& deltaH_Table,
                              Zuzax::DenseMatrix& deltaS_Table,
                              Zuzax::DenseMatrix& Afwd_Table,
                              Zuzax::DenseMatrix& EafwddivR_Table,
                              Zuzax::DenseMatrix& Arev_Table,
                              Zuzax::DenseMatrix& EarevdivR_Table,
                              Zuzax::DenseMatrix& kfwdPrime_Table,
                              Zuzax::DenseMatrix& krevPrime_Table);

void getGERKineticsTables(TemperatureTable& TT, Zuzax::PhaseList* pl, Zuzax::Kinetics& kin, Zuzax::ExtraGlobalRxn& egr,
                     RxnTempTableStuff& rts);

void printGERKineticsTable(Zuzax::PhaseList* pl, int j, TemperatureTable& TT,
                           Zuzax::Kinetics& kin, Zuzax::ExtraGlobalRxn& egr, Zuzax::RxnMolChange* rmc,
                           RxnTempTableStuff& rts);

void processGERCurrentVsPotTable(Zuzax::RxnMolChange* rmc, Zuzax::PhaseList* pl, int irxn,
                                 TemperatureTable& TT,
                                 Zuzax::Kinetics& kin,
                                 Zuzax::ExtraGlobalRxn& egr,
                                 RxnTempTableStuff& rts);

void printAffinityHeader(Zuzax::RxnMolChange* rmc,
                         Zuzax::PhaseList* pl,
                         int irxn,
                         TemperatureTable& TT,
                         Zuzax::Kinetics& kin,
                         Zuzax::DenseMatrix& kfwd_Table,
                         Zuzax::DenseMatrix& krev_Table,
                         Zuzax::DenseMatrix& deltaG_Table,
                         Zuzax::DenseMatrix& deltaH_Table,
                         Zuzax::DenseMatrix& deltaS_Table,
                         Zuzax::DenseMatrix& Afwd_Table,
                         Zuzax::DenseMatrix& EafwddivR_Table,
                         Zuzax::DenseMatrix& Arev_Table,
                         Zuzax::DenseMatrix& EarevdivR_Table,
                         Zuzax::DenseMatrix& kfwdPrime_Table,
                         Zuzax::DenseMatrix& krevPrime_Table,
                         double* unitsfwd,
                         double* unitsRev);

void processAffinityTable(Zuzax::RxnMolChange* rmc,
                          Zuzax::PhaseList* pl,
                          int irxn,
                          TemperatureTable& TT,
                          Zuzax::Kinetics& kin,
                          Zuzax::DenseMatrix& kfwd_Table,
                          Zuzax::DenseMatrix& krev_Table,
                          Zuzax::DenseMatrix& deltaG_Table,
                          Zuzax::DenseMatrix& deltaH_Table,
                          Zuzax::DenseMatrix& deltaS_Table,
                          Zuzax::DenseMatrix& Afwd_Table,
                          Zuzax::DenseMatrix& EafwddivR_Table,
                          Zuzax::DenseMatrix& Arev_Table,
                          Zuzax::DenseMatrix& EarevdivR_Table,
                          Zuzax::DenseMatrix& kfwdPrime_Table,
                          Zuzax::DenseMatrix& krevPrime_Table);

//==================================================================================================================================

#endif
