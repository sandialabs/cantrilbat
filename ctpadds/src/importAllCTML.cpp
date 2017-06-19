/**
 *  @file importAllCTML.cpp
 *     This file contains routines which are global routines, i.e., not part of any object,which take as input, ctml
 *     pointers to data, and pointers to Zuzax objects in order to initialize the Zuzax objects with data 
 *     from the ctml tree structures.
 */

/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#include "cantera/base/ctml.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/StoichSubstanceSSTP.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"

#include "importAllCTML.h"
#include "ElectrolyteKinetics.h"
#include "SolidKinetics.h"

#include <cstring>

#ifdef useZuzaxNamespace
using namespace ztml;
#else
using namespace ctml;
#endif

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
thermo_t_double* processExpandedThermoPhase(XML_Node* const xmlphase)
{
    thermo_t_double* tPhase = nullptr;
    try {
        tPhase = Zuzax::newPhase(*xmlphase);
    } catch (UnknownThermoPhaseModel& uName) {
        const XML_Node& th = xmlphase->child("thermo");
        std::string model = th["model"];
        printf("\t\t model = %s\n", model.c_str());
        /*
         *  IdealSolidSolnPhase: This is not yet in the Factory, but is in the main Zuzax distribution.
         */
        if (model == "IdealSolidSolution") {
            popError();
            tPhase = new IdealSolidSolnPhase(*xmlphase);
        }
        /*
         * IdealMolalSoln: This is not yet in the Factory, but is in the main Zuzax distribution.
         */
        else if (model == "IdealMolalSoln") {
            popError();
            tPhase = new IdealMolalSoln(*xmlphase);
        }
        /*
         * DebyeHuckel:    This is not yet in the Factory, but is in the main Zuzax distribution.
         */
        else if (model == "DebyeHuckel") {
            popError();
            tPhase = new DebyeHuckel(*xmlphase);
        }
        /*
         * HMWSoln:        This is not yet in the Factory, but is in the main Zuzax distribution.
         */
        else if (model == "HMW" || model == "HMWSoln") {
            popError();
            tPhase = new HMWSoln(*xmlphase);
        }
        /*
         * StoichSubstance: This is not yet in the Factory, but is in the main Zuzax distribution.
         */
        else if (model == "StoichSubstance") {
            popError();
            tPhase = new StoichSubstanceSSTP(*xmlphase);
        } else {
            throw;
        }
    }
    return tPhase;
}
//==================================================================================================================================
Kinetics* processExpandedKinetics(XML_Node* const xmlPhase, std::vector<thermo_t_double*> tpList)
{
    Kinetics* kin = nullptr;
    if (tpList.size() == 0) {
        throw ZuzaxError("processExpandedKinetics()", "Expecting tpList to be filled before entry\n");
    }
    thermo_t_double* tp = tpList[0];
    if (!tp) {
        throw ZuzaxError("processExpandedKinetics()", "Expecting tpList to be filled before entry\n");
    }

    XML_Node& kinNode = xmlPhase->child("kinetics");
    std::string kModel = kinNode.attrib("model");
    if (kModel == "" || (!strcasecmp(kModel.c_str(), "None"))) {
        return nullptr;
    }
    /*
     * SolidKinetics:       This is a New one.  It's for solid phases
     */
    else if (kModel == "SolidKinetics") {
        SolidKinetics* sk_ptr = new SolidKinetics(tp);
        sk_ptr->importMechanism(*xmlPhase, "");
        kin = sk_ptr;
    } else {
        try {
            kin = newKineticsMgr(*xmlPhase, tpList);
        } catch (UnknownKineticsModel) {
            XML_Node& kinNode = xmlPhase->child("kinetics");
            std::string kModel = kinNode.attrib("model");
            if (kModel == "SolidKinetics") {
                SolidKinetics* sk_ptr = new SolidKinetics(tp);
                sk_ptr->importMechanism(*xmlPhase, "");
                kin = sk_ptr;
                popError();
            } else if (kModel == "Electrolyte") {
                ElectrolyteKinetics* ek_ptr = new ElectrolyteKinetics(tp);
                ek_ptr->importMechanism(*xmlPhase, "");
                kin = ek_ptr;
            } else if (kModel == "NONE") {
                popError();
            } else {
                showErrors(std::cout);
            }
        }
    }
    return kin;
}
//==================================================================================================================================
InterfaceKinetics* processExpandedInterfaceKinetics(XML_Node* const xmlPhase, std::vector<thermo_t_double*> tpList)
{
    InterfaceKinetics* kin = nullptr;
    if (tpList.size() == 0) {
        throw ZuzaxError("processExpandedKinetics()", "Expecting tpList to be filled before entry");
    }
    thermo_t_double* tp = tpList[0];
    if (!tp) {
        throw ZuzaxError("processExpandedKinetics()", "Expecting tpList to be filled before entry\n");
    }
    XML_Node& kinNode = xmlPhase->child("kinetics");
    std::string kModel = kinNode.attrib("model");
    if (kModel == "" || (!strcasecmp(kModel.c_str(), "None"))) {
        return nullptr;
    } else {
        try {
            Kinetics* kinBase = newKineticsMgr(*xmlPhase, tpList);
            kin = dynamic_cast<InterfaceKinetics*>(kinBase);
            if (kin == 0) {
                throw ZuzaxError("ProcessExpandedInterfaceKinetics()", " Dynamic cast to InterfaceKinetics failed");
            }
        } catch (UnknownKineticsModel& ee) {
            XML_Node& kinNode = xmlPhase->child("kinetics");
            std::string kModel = kinNode.attrib("model");
            if (kModel == "Electrode") {
                popError();
                throw ZuzaxError(" processExpandedInterfaceKinetics", "Electrode model selected");
                kin = 0;
                popError();
            } else if (kModel == "NONE") {
                kin = 0;
                popError();
            } else {
                showErrors(std::cout);
                throw ee;
            }
        }
    }
    return kin;
}
//==================================================================================================================================
Transport* processExpandedTransport(const XML_Node* const xmlPhase, thermo_t_double* const tp)
{
    Transport* tran = nullptr;
    if (xmlPhase->hasChild("transport")) {
        XML_Node& tranNode = xmlPhase->child("transport");
        std::string tModel = tranNode.attrib("model");
        if (tModel == "" || (!strcasecmp(tModel.c_str(), "None"))) {
            return tran;
        } else {
            tran = newTransportMgr(tModel, tp);
        }
    }
    return tran;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

