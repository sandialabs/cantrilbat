/*
 *
 *  xmlSolnDiff File1.xml File2.xml
 *
 *  Compares the variable values in two Zuzax solution xml
 *  files.
 *  The comparison is done using a weighted norm basis.
 *
 *  The two files should be basically equal. However, File1.xml is
 *  taken as the reference file, that has precedence, when there is
 *  something to be decided upon.
 *
 *  Arguments:
 *   -h = prints this usage information
 *
 *  Shell Return Values
 *    1 = Comparison was successful
 *    0 = One or more nodal values failed the comparison
 *   -1 = Apples to oranges, the files can not even be compared against
 *        one another.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <fstream>
#include <unistd.h>


#include "zuzax/base/xml.h"
#include "zuzax/base/ctml.h"
#include "zuzax/base/ctexceptions.h"
#include "zuzax/base/zzcompare.h"
#include "zuzax/base/VarType.h"

#include "zuzax/numerics/SolnOutput.h"
#include "zuzax/numerics/SolnState.h"
#include "zuzax/numerics/TimeSpec.h"
#include "zuzax/numerics/GridSpec.h"
#include "zuzax/base/Distrib1D.h"

#include "zuzax/kinetics/Kinetics.h"
#include "zuzax/kinetics/ReactionData.h"
#include "zuzax/kinetics/importKinetics.h"

using namespace std;
using namespace Zuzax;
using namespace ztml;
//==================================================================================================================================
//! Read an XML file into a XML_Node Tree structure
/*!
 *  @param[in]               inputFile           Name of the input file
 */
static XML_Node* readXML(std::string inputFile)
{

    if (inputFile.size() == 0) {
        throw ZuzaxError("readXML()",  "input file is null");
    }
    std::string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
        throw ZuzaxError("readXML()","could not open " +path+" for reading.");
    }
    XML_Node* fxml = new XML_Node();
    fxml->build(fin);
    return fxml;
}

//==================================================================================================================================
static XML_Node* getReactionDataXML(XML_Node* xmlTop, int num = 0)
{
    XML_Node* xctml = xmlTop;
    if (xctml->name() != "ctml") {
        xctml = xmlTop->findByName("ctml");
        if (!xctml) {
            throw ZuzaxError("getSimulNum","Can't find ctml node");
        }
    }

    XML_Node* pp;
    if (xctml->hasChild("reactionData")) {
         pp = & xctml->child("reactionData");
    } else {
        throw ZuzaxError("readXML()","XML file doesn't have a phase in it");
    }
    return pp;
}
//==================================================================================================================================
static void getStick(const XML_Node& node, ReactionData& r, doublevalue& A, doublevalue& b, doublevalue& E)
{
    /*
     * species is the name of the special reactant whose surface flux rate will be calculated.
     *      isp = species # in the local phase
     *      ispKinetics = species # in the kinetics object
     *      ispPhaseIndex = phase # of the special species
     */
    std::string spname = node["species"];

    // Gather the molecular weight of the sticking coefficient species

    A = getFloat(node, "A", "toSI");
    b = getFloat(node, "b");
    E = getFloat(node, "E", "actEnergy");

}
//==================================================================================================================================
static void getArrhenius(const XML_Node& node, int& labeled, doublevalue& A, doublevalue& b, doublevalue& E, doublevalue& T0)
{
    if (node["name"] == "k0") {
        labeled = -1;
    } else if (node["name"] == "kHigh") {
        labeled = 1;
    } else {
        labeled = 0;
    }
    /*
     * We parse the children for the A, b, and E components. And we parse for the alternative k0 form adding a T0.
     */
    T0 = 0.0;
    if (node.hasChild("A")) {
        A = ztml::getFloat(node, "A", "toSI");
        b = ztml::getFloat(node, "b");
        E = ztml::getFloat(node, "E", "actEnergy");
    } else if (node.hasChild("k0")) {
        A = ztml::getFloat(node, "k0", "toSI");
        b = ztml::getFloat(node, "b");
        E = ztml::getFloat(node, "E", "actEnergy");
        T0 = ztml::getFloat(node, "T0");
        if (T0 <= 0.0) {
            throw ZuzaxError("importKinetics::getArrhenius()", "T0 must be positive. It's a temperature");
        }
        // Convert k0 to a regular preexponential A

        A *= exp(E/T0) / pow(T0,b) ;
    } else {
        throw ZuzaxError("importKinetics::getArrhenius()", "Arrhenius node doesn't have either an A or an k0 element");
    }
}
//==================================================================================================================================
static int getRxnsXML(XML_Node* xmlRD, std::vector<XML_Node*>& ccc)
{
    xmlRD->getChildren("reaction", ccc);
    int sz = ccc.size();
    return sz;
}
//==================================================================================================================================
static bool getReagentsB(ReactionData& rdata, const XML_Node& rxnNode, const ReactionRules& rules)
{
    /*
     * The id of reactants and products are kept in child elements
     * of reaction, named "reactants" and "products". 
     */
    const XML_Node& rgr = rxnNode.child("reactants");
    const XML_Node& rgp = rxnNode.child("products");
    rdata.reactantsCompMap = parseCompString(rgr.value());
    rdata.productsCompMap = parseCompString(rgp.value());
    return true;
}
//==================================================================================================================================
static bool getRateCoefficient(const XML_Node& kf, ReactionData& rdata, const ReactionRules& rules)
{
    if (rdata.reactionType == PLOG_RXN) {
        rdata.rateCoeffType = PLOG_REACTION_RATECOEFF_TYPE;
        for (size_t m = 0; m < kf.nChildren(); m++) {
            const XML_Node& node = kf.child(m);
            doublevalue p = ztml::getFloat(node, "P", "toSI");
            vector_fp& rate = rdata.plogParameters.insert( std::make_pair(p, vector_fp()))->second;
            rate.resize(3);
            rate[0] = ztml::getFloat(node, "A", "toSI");
            rate[1] = ztml::getFloat(node, "b");
            rate[2] = ztml::getFloat(node, "E", "actEnergy") / GasConstant;
        }

    } else if (rdata.reactionType == CHEBYSHEV_RXN) {
        rdata.rateCoeffType = CHEBYSHEV_REACTION_RATECOEFF_TYPE;
        rdata.chebTmin = getFloat(kf, "Tmin", "toSI");
        rdata.chebTmax = getFloat(kf, "Tmax", "toSI");
        rdata.chebPmin = getFloat(kf, "Pmin", "toSI");
        rdata.chebPmax = getFloat(kf, "Pmax", "toSI");
        const XML_Node& coeffs = kf.child("floatArray");
        rdata.chebDegreeP = atoi(coeffs["degreeP"].c_str());
        rdata.chebDegreeT = atoi(coeffs["degreeT"].c_str());
        getFloatArray(kf, rdata.chebCoeffs, false);

    } else {

        std::string type = kf.attrib("type");
        if (type == "") {
            type = "Arrhenius";
            rdata.rateCoeffType = ARRHENIUS_REACTION_RATECOEFF_TYPE;
        }
        if (type == "ExchangeCurrentDensity") {
            rdata.rateCoeffType = EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE;
        } else if (type == "Arrhenius") {

        } else {
            throw ZuzaxError("getRateCoefficient", "Unknown type: " + type);
        }

        vector_fp c_alt(4,0.0), c_base(4,0.0);
        for (size_t m = 0; m < kf.nChildren(); m++) {
            const XML_Node& c = kf.child(m);
            std::string nm = c.name();
            int labeled=0;

            if (nm == "Arrhenius") {
                vector_fp coeff(4);
                if (c["type"] == "stick") {
                    getStick(c, rdata, coeff[0], coeff[1], coeff[2]);
                    c_base = coeff;
                } else {
                    getArrhenius(c, labeled, coeff[0], coeff[1], coeff[2], coeff[3]);
                    if (labeled == 0 || rdata.reactionType == THREE_BODY_RXN
                            || rdata.reactionType == ELEMENTARY_RXN) {
                        c_base = coeff;
                    } else {
                        c_alt = coeff;
                    }
                }
                if (rdata.reactionType == SURFACE_RXN || rdata.reactionType == EDGE_RXN) {
                    //getCoverageDependence(c, kin.thermo(kin.reactionPhaseIndex()), rdata);
                    throw ZuzaxError("  get Rate", "getConverageDep not covered");
                 }

                if (coeff[0] < 0.0 && !rules.allowNegativeA) {
                    throw ZuzaxError("getRateCoefficient",
                                       "negative A coefficient for reaction "+int2str(rdata.number));
                }


           } else if (nm == "Arrhenius_ExchangeCurrentDensity") {
                vector_fp coeff(4);
                getArrhenius(c, labeled, coeff[0], coeff[1], coeff[2], coeff[3]);
                c_base = coeff;
                rdata.rateCoeffType = EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE;
            } else if (nm == "falloff") {
                //getFalloff(c, rdata);
                throw ZuzaxError("getRates", "Falloff reactions not covered");
            } else if (nm == "efficiencies") {
                //getEfficiencies(c, kin, rdata, rules);
                 throw ZuzaxError("getRates", "getEfficiencies not covered");
            } else if (nm == "electrochem") {
                rdata.beta = fpValue(c["beta"]);
            }
        }
        /*
         * Store the coefficients in the ReactionData object for return
         * from this function.
         */
        if (rdata.reactionType == FALLOFF_RXN) {
            rdata.rateCoeffParameters = c_base;
            rdata.auxRateCoeffParameters = c_alt;
        } else if (rdata.reactionType == CHEMACT_RXN) {
            rdata.rateCoeffParameters = c_alt;
            rdata.auxRateCoeffParameters = c_base;
        } else {
            rdata.rateCoeffParameters = c_base;
        }

    }
    return true;
}
//==================================================================================================================================
//! We read enough of the reaction data information to print out what we need in the table.
//! No checking of anything is done.
static int readReactionDataXML(int iRxn, const XML_Node& rxnNode,  ReactionData& rdata)
{

    ReactionRules rules;
      /*
     * Search the reaction element for the attribute "type".
     * If found, then branch on the type, to fill in appropriate
     * fields in rdata.
     */
    rdata.reactionType = ELEMENTARY_RXN;
    std::string typ = rxnNode["type"];
    std::string ltype = lowercase(typ);
    if (typ == "falloff") {
        rdata.reactionType = FALLOFF_RXN;
        rdata.falloffType = SIMPLE_FALLOFF;
    } else if (typ == "chemAct") {
        rdata.reactionType = CHEMACT_RXN;
        rdata.falloffType = SIMPLE_FALLOFF;
    } else if (typ == "threeBody") {
        rdata.reactionType = THREE_BODY_RXN;
    } else if (typ == "plog") {
        rdata.reactionType = PLOG_RXN;
    } else if (typ == "chebyshev") {
        rdata.reactionType = CHEBYSHEV_RXN;
    } else if (typ == "surface") {
        rdata.reactionType = SURFACE_RXN;
    } else if (typ == "edge") {
        rdata.reactionType = EDGE_RXN;
    } else if (ltype == "butlervolmer_noactivitycoeffs") {
        rdata.reactionType = BUTLERVOLMER_NOACTIVITYCOEFFS_RXN;
    } else if (ltype == "butlervolmer_constantcurrentdensity") {
        rdata.reactionType = BUTLERVOLMER_CONSTANTCURRENTDENSITY_RXN;
    } else if (ltype == "butlervolmer") {
        rdata.reactionType = BUTLERVOLMER_RXN;
    } else if (ltype == "surfaceaffinity") {
        rdata.reactionType = SURFACEAFFINITY_RXN;
    } else if (ltype == "global") {
        rdata.reactionType = GLOBAL_RXN;
    } else if (typ != "") {
        throw ZuzaxError("installReaction()", "Unknown reaction type: " + typ);
    }

  // Read the id field. If n
    std::string idRxn = rxnNode["id"];
    if (idRxn == "") {
        rdata.RxnID = fp2str(iRxn);
    } else {
        rdata.RxnID = idRxn;
    }

    // Check to see if the reaction is specified to be a duplicate of another
    // reaction. It's an error if the reaction is a duplicate and this is not set.
    rdata.duplicate = (rxnNode.hasAttrib("duplicate")) ? 1 : 0;

    // Check to see if the reaction rate constant can be negative. It's an
    // error if a negative rate constant is found and this is not set.
    rules.allowNegativeA = (rxnNode.hasAttrib("negative_A")) ? 1 : 0;

    // Use the contents of the "equation" child element as the reaction's
    // string representation. Post-process to convert "[" and "]" characters
    // back into "<" and ">" which cannot easily be stored in an XML file. This
    // reaction string is used only for display purposes. It is not parsed for
    // the identities of reactants or products.
    rdata.equation = (rxnNode.hasChild("equation")) ? rxnNode("equation") : "<no equation>";
    static const char* delimiters[] = {" [=] ", " =] ", " = ", "[=]", "=]", "="};
    static const char* replacements[] = {" <=> ", " => ", " = ", "<=>", "=>", "="};
    for (size_t i = 0; i < 6; i++) {
        size_t n = rdata.equation.find(delimiters[i]);
        if (n != npos) {
            size_t w = strlen(delimiters[i]);
            rdata.reactantString = stripws(rdata.equation.substr(0, n));
            rdata.productString = stripws(rdata.equation.substr(n+w, npos));
            rdata.equation.replace(n, w, replacements[i]);
            break;
        }
    }


     getReagentsB(rdata, rxnNode, rules);

    XML_Node& rateXML =  rxnNode.child("rateCoeff");

    getRateCoefficient(rateXML, rdata, rules);


    return 0;
}
//==================================================================================================================================
int main(int argc, char* argv[])
{
    char buf[100000];
    const char*  fileName1=NULL;
    XML_Node* xmlTop = nullptr;
    std::string optArgSS, arglc;
    std::vector<XML_Node*> ccc;
    std::string src;
    fileName1 = "aqueousRxns.xml"; 
    /*
     * Interpret command line arguments
     */
    if (!(xmlTop = readXML(fileName1))) {
        fprintf(stderr,"Error opening up file1, %s\n", fileName1);
        exit(-1);
    }


    Zuzax::XML_Node* xmlRD = getReactionDataXML(xmlTop);

    size_t numParamBlocks = getRxnsXML(xmlRD, ccc);

    FILE *fcsv = fopen("RxnTable.csv", "w");    
    printf("Got  %d rxn points\n", (int) numParamBlocks);

    std::vector<ReactionData> rdataVec(numParamBlocks);

    // Buffer the csv line into buf;
    for (size_t iP = 0; iP < numParamBlocks; iP++) {
        XML_Node* rrNode = ccc[iP];
        buf[0] = '\0';
        int slen = 0;
        ReactionData& rdata = rdataVec[iP];

        int notOK =  readReactionDataXML(iP, *rrNode,  rdata);
        if (notOK) {
            throw ZuzaxError("readData", "Rxn data can't be handled by this limited interpretor");
        }


        sprintf(buf + slen, " %s , " , rdata.equation.c_str());
        slen = strlen(buf);

        sprintf(buf + slen, " %g , %g , %g , ", rdata.rateCoeffParameters[0],  rdata.rateCoeffParameters[1],
                rdata.rateCoeffParameters[2]);
        slen = strlen(buf);

        const XML_Node* xsource = rrNode->findByName("source");
        src = "";
        if (xsource) {
            src = xsource->value();
        }
        sprintf(buf + slen, " \"%s\" ",  src.c_str());
        slen = strlen(buf);
 
        fprintf(fcsv,"%s\n", buf);

    }

    fclose (fcsv);

    return 0;
}
//==================================================================================================================================
