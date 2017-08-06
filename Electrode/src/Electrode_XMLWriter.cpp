/**
 *  @file Electrode_XMLWriter.cpp
 *     Implementation of the XML read and write routines.
 *     (see \ref electrode_mgr and class \link Zuzax::Electrode Electrode\endlink).
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "Electrode.h"
#include <fstream>

using namespace std;
//-----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//===================================================================================================================================
// Specifies the amount of output that the Electrode object writes to its solution file
/*
 *    The level is given by the following table
 *            - 0     Zero output (default)
 *            - 1     Global time end
 *            - 2     global time increment
 *            - 3     intermediate steps
 */
void Electrode::specifySolutionFileLevel(int level, const char* const baseName)
{
    printXMLLvl_ = level;
    baseNameSoln_ = baseName;
}
//===================================================================================================================================
//  Create a timeIncrement XML element to store the results for intermediate steps of the solver.
/*
 *    Creates the following XML tree structure.
 *    If the boolean addInitState is false the timeState t_init is not written.
 *
 *     <timeIncrement   type="intermediate">
 *        <timeState type="t_init">
 *              ....xmlStateData_init_
 *        </timeState>
 *        <timeState type="t_intermediate">
 *              ....xmlStateData_final_
 *        </timeState>
 *     </timeIncrement>
 *
 *    This XML Tree is storred in the variable  xmlTimeIncrementIntermediateData_
 *
 *  @param addInitState   Boolean that if true adds the initial state to the tree.
 *                        The default is true.
 */
void Electrode::makeXML_TI_intermediate(bool addInitState)
{
    std::string fmt = "%22.14E";
    delete xmlTimeIncrementIntermediateData_;
    xmlTimeIncrementIntermediateData_ = new XML_Node("timeIncrement");
    xmlTimeIncrementIntermediateData_->addAttribute("type", "intermediate");
    if (addInitState) {
        XML_Node* xmi = new XML_Node("timeState");
        xmi->addAttribute("type", "t_init");
        xmi->addAttribute("domain", electrodeDomainNumber_);
        xmi->addAttribute("cellNumber", electrodeCellNumber_);
        xmi->addChild("time", tinit_, fmt);
        xmi->addChild(*xmlStateData_init_);
        xmlTimeIncrementIntermediateData_->mergeAsChild(*xmi);
    }
    XML_Node* xmt = new XML_Node("timeState");
    xmt->addAttribute("type", "t_intermediate");
    xmt->addAttribute("domain", electrodeDomainNumber_);
    xmt->addAttribute("cellNumber", electrodeCellNumber_);
    xmt->addChild("time", tfinal_, fmt);
    xmt->addChild(*xmlStateData_final_);
    xmlTimeIncrementIntermediateData_->mergeAsChild(*xmt);
}

//==================================================================================================================================
// Creates a timeIncrement XML element to store the results for global steps of the Electrode solver
/*
 *   Creates a XML Tree structure of the following form. What this function does is to delete the old record
 *   and start a new record. Later calls to addtoXML_TI_final() adds records to the XML tree.
 *
 *     <timeIncrement   type="global">
 *        <timeState type="t_init">
 *              ....xmlStateData_init_
 *        </timeState>
 *        <timeState type="t_intermediate">
 *              ....xmlStateData_final_
 *        </timeState>
 *        <timeState type="t_inal">
 *              ....xmlStateData_final_
 *        </timeState>
 *     </timeIncrement>
 *
 *    This XML Tree is storred in the variable  xmlTimeIncrementData_
 *
 *
 *  @param addInitState   Boolean that if true adds the initial state to the tree.
 *                        The default is true.
 */
void Electrode::startXML_TI_final(bool addInitState)
{
    std::string fmt = "%22.14E";
    static int firstTime = 0;
    delete xmlTimeIncrementData_;
    xmlTimeIncrementData_ = new XML_Node("timeIncrement");
    xmlTimeIncrementData_->addAttribute("type", "global");
    if (firstTime > 0 && printXMLLvl_ <= 1) {
        addInitState = false;
    }
    if (addInitState) {
        XML_Node* xmi = new XML_Node("timeState");
        xmi->addAttribute("type", "t_init");
        xmi->addAttribute("domain", electrodeDomainNumber_);
        xmi->addAttribute("cellNumber", electrodeCellNumber_);
        xmi->addChild("time", tinit_, fmt);
        xmi->addChild(*xmlStateData_init_);
        xmlTimeIncrementData_->mergeAsChild(*xmi);
    }
}
//===================================================================================================================================
// Adds to a timeIncrement XML element to store the results for intermediate or global-final steps of the solver.
/*
 *    Adds a timeState record to the following XML tree structure.
 *
 *     <timeIncrement   type="global">
 *        <timeState type="t_init">
 *              ....xmlStateData_init_
 *        </timeState>
 *        <timeState type="t_intermediate">
 *              ....xmlStateData_final_
 *        </timeState>
 *        <timeState type="t_final">
 *              ....xmlStateData_final_
 *        </timeState>
 *     </timeIncrement>
 *
 *    This XML Tree is storred in the variable  xmlTimeIncrementData_
 *
 *  @param notDone   Boolean that if true sets the type attribute to t_intermediate
 *                   If false, the type attribute is set to t_final
 */
void Electrode::addtoXML_TI_final(bool notDone)
{
    std::string fmt = "%22.14E";
    // leave out t_intermediate records for some print levels
    if (notDone) {
        if (printXMLLvl_ == 1 || printXMLLvl_ == 2) {
            return;
        }
    }
    XML_Node* xmt = new XML_Node("timeState");
    if (notDone) {
        xmt->addAttribute("type", "t_intermediate");
    } else {
        xmt->addAttribute("type", "t_final");
    }
    xmt->addAttribute("domain", electrodeDomainNumber_);
    xmt->addAttribute("cellNumber", electrodeCellNumber_);
    xmt->addChild("time", tfinal_, fmt);
    xmt->addChild(*xmlStateData_final_);

    xmlTimeIncrementData_->mergeAsChild(*xmt);
}
//===================================================================================================================================
//  Returns true if the step was successful and false otherwise
bool Electrode::writeTimeStateFinal_toXML(XML_Node& bb)
{
    std::string fmt = "%22.14E";
    
    // There may be situations when we haven't created the state data. When that happens we will create the
    // the state data right here. This currently occurs for writing out initial conditions.
    if (!xmlStateData_final_) {
	if (eState_save_) {
	    xmlStateData_final_ = eState_save_->write_electrodeState_ToXML();
	}
    }
    if (xmlStateData_final_) {
	XML_Node&ts = bb.addChild("timeState");
	ts.addAttribute("type", "t_final");
	ts.addAttribute("domain", electrodeDomainNumber_);
	ts.addAttribute("cellNumber", electrodeCellNumber_);
	ts.addChild("time", tfinal_, fmt);
	
	//  Add information from the saved solution to the record. (Lengthy operation)
	ts.addChild(*xmlStateData_final_);
	return true;
    }
    return false;
}
//==================================================================================================================================
// Given a Time increment record this routine loads the saved solution for t_final into the electrode object
/*
 *
 *     This routine reads an XML tree layout, such as the example below. It is given the globalTimeStep XML element
 *     pointer. It then reads the timeState  final record and final time, and initializes the current object with
 *     the state contained in the record.
 *
 *            <globalTimeStep index = 312>
 *               <timeIncrement type="global">
 *                  <timeState type="t_init">
 *                     ....xmlStateData_init_
 *                  </timeState>
 *                  <timeState type="t_intermediate">
 *                      ....xmlStateData_final_
 *                  </timeState>
 *                  <timeState type="t_final">
 *                     <time> 3.45 <time>               <-----------  time used.
 *                     <electrodeState>                 <-----------  Record read
 *                     </electrodeState>
 *                  </timeState>
 *               </timeIncrement>
 *            </globalTimeStep>
 *
 *
 *
 *  @param xGSTI    Global time step increment record
 */
void Electrode::loadGlobalTimeStepTFinalState(const XML_Node* const xGTS)
{
    // Check name
    std::string ss = xGTS->name();
    if (ss != "globalTimeStep") {
        throw Electrode_Error("Electrode::loadGlobalTimeStepTFinalState()", "Expected node named globalTimeStep. got " + ss);
    }
    const XML_Node* const xGTSI = xGTS->findByName("timeIncrement");
    if (!xGTSI) {
        throw Electrode_Error("Electrode::loadGlobalTimeStepTFinalState()", "could not find a node named timeIncrement");
    }
    if (pendingIntegratedStep_) {
        throw Electrode_Error("Electrode::setTime()", "called when there is a pending step");
    }
    /*
     *   Search the immediate children for the t_final attribute.
     *   Store the timeState XML record's pointer in rTFinal
     */
    const XML_Node* const rTFinal = xGTSI->findByAttr("type", "t_final", 1);
    if (!rTFinal) {
        throw Electrode_Error("Electrode::loadGlobalTimeStepTFinalState()",
                              "could not find a node with attribute type t_final within timeIncrement");
    }
    // load
    loadTimeState(*rTFinal);
}
//==================================================================================================================================
void Electrode::loadGlobalTimeStepTInitState(const XML_Node* const xGTS)
{
    // Check name
    std::string ss = xGTS->name();
    if (ss != "globalTimeStep") {
        throw Electrode_Error("Electrode::loadGlobalTimeStepTInitState()", "Expected node named globalTimeStep. got " + ss);
    }
    const XML_Node* const xGTSI = xGTS->findByName("timeIncrement");
    if (!xGTSI) {
        throw Electrode_Error("Electrode::loadGlobalTimeStepTInitState()", "could not find a node named timeIncrement");
    }
    if (pendingIntegratedStep_) {
        throw Electrode_Error("Electrode::loadGlobalTimeStepTInitState()", "called when there is a pending step");
    }
    
    //   Search the immediate children for the t_init attribute.
    //   Store the timeState XML record's pointer in rTInit
    const XML_Node* const rTInit = xGTSI->findByAttr("type", "t_init", 1);
    if (!rTInit) {
        throw Electrode_Error("Electrode::loadGlobalTimeStepTInitState()",
                              "could not find a node with attribute type t_init within timeIncrement");
    }
    
    loadTimeState(*rTInit);
}
//==================================================================================================================================
double Electrode::loadTimeState(const XML_Node& xTimeState)
{
    //  Check name of the XML element that will be read to make sure we're in the right location
    std::string ss = xTimeState.name();
    if (ss != "timeState") {
        throw Electrode_Error("Electrode::loadTimeState()", "Expecting the XML element \"timeState\", but got %s", ss.c_str());
    }
    std::string type = xTimeState["type"];
    
    //  Get the time of the solution saved in the 
    double time = ZZctml::getFloat(xTimeState, "time");
    
    //  Get the XML electrodeState record from the timeState XML element.
    const XML_Node* const xState = xTimeState.findByName("electrodeState");
    if (!xState) {
        throw Electrode_Error("Electrode::loadTimeState()", "Could not find the electrodeState XMl element within XML state record");
    }

    //  Read the state of the electrode into the Estate object
    eState_save_->readStateFromXML(*xState);

    //  Then, set this Electrode's internal state from eState current values. 
    //  When we do this all init and final states are set to this state.
    eState_save_->setStateElectrode_fromEState(this);
    
    //  Set the time within the electrode object
    setTime(time);
    return time;
}
//==================================================================================================================================
//  Select the global time step increment record by the consequatively numbered record index number
XML_Node* Electrode::selectGlobalTimeStepIncrement(XML_Node* xSoln, int globalTimeStepNum)
{
    //  Find the electrodeOutput XML element.
    //     -> Later we will generalize this to search amongst multiple electrodeOutput objects
    XML_Node* eOutput = xSoln->findByName("electrodeOutput");
    if (!eOutput) {
        throw Electrode_Error("Electrode::selectGlobalTimeStepIncrement()", "Could not find a node named electrodeOutput");
    }
    //  Search for a particular global step number
    XML_Node* eRecord = eOutput->findNameIDIndex("globalTimeStep", "", globalTimeStepNum);
    if (!eRecord) {
        throw Electrode_Error("Electrode::selectGlobalTimeStepIncrement()", "Could not find a node named globalTimeStep");
    }
  
    //  Return a pointer to the record
    return eRecord;
}
//==================================================================================================================================
//  Write a global time step to the solution file. 
//  Wrap the timeIncrement XML element within a solution XML element and then write it out to an output file
/*
 *   We assume that the XML tree has the following topology
 *
 *         <ctml>
 *           <electrodeOutput index = 1>
 *             <timeStamp>
 *                 Jan 19, 2012
 *             </timeStamp>
 *             <globalTimeStep index = 1>
 *               <timeIncrement type="global">
 *                 ...
 *               </timeIncrement>
 *             </globalTimeStep>
 *             <globalTimeStep index = 2>
 *               <timeIncrement type="global">
 *                 ...
 *               </timeIncrement>
 *             </globalTimeStep>
 *             . . .
 *           </electrodeOutput>
 *         </ctml>
 *
 *  This routine adds another <timeIncrement> XML element to the end of the file. It
 *  takes care to first eliminate any existing to backspace over the last
 *  </electrodeOutput> and </ctml> entries before writing the new <timeIncrement> XML  element.
 */
void Electrode::writeSolutionTimeIncrement(bool startNewRecord, bool reset, int stepNumOverride)
{
    static bool haventBeenHere = true;
    bool doTrunc = false;

    //  solnNum is the global solution number. There can be more than one time dependent solution in an output file.
    static int solnNum = 1;

    if (!eState_save_) {
        throw Electrode_Error("Electrode::writeSolutionTimeIncrement()", 
                              "Solution Output has been requested, but it has never been set up by invoking electrode_stateSave_create()");
    }
    /*
     *  stepNum is the global time step number written to the file.
     */
    static int stepNum = 0;

    if (haventBeenHere) {
        solnNum = 1;
        stepNum = 0;
        doTrunc = true;
    } else if (startNewRecord) {
        solnNum++;
        stepNum = 0;
    }
    if (reset) {
        haventBeenHere = true; 
        solnNum = 1;
        stepNum = 0;
        doTrunc = true;
    }

    time_t aclock;
    ::time(&aclock); /* Get time in seconds */
    struct tm*   newtime = localtime(&aclock); /* Convert time to struct tm form */

    /*
     *  Determine the name of the solution file
     */
    std::string bname = "soln";
    if (baseNameSoln_ != "") {
        bname = baseNameSoln_;
    }
    /*
     *  Add in the domain number and cell number
     */
    std::string fname = (bname + "_" + int2str(electrodeDomainNumber_) + "_" + int2str(electrodeCellNumber_) + ".xml");
    /*
     *  Start a new root XML element and put the top elements in
     */
    ZZCantera::XML_Node root("--");
    ZZCantera::XML_Node& ct = root.addChild("ctml");
    ZZCantera::XML_Node& soln = ct.addChild("electrodeOutput");

    /*
     *  Add the integer attribute, index, to the electrodeOutput element.
     *  -> so far we have an index of one. However, if there are continuation runs in the
     *     file and we want to add more than one
     */
    std::string ii = ZZCantera::int2str(solnNum);
    soln.addAttribute("index", ii);

    stepNum++;
    if (stepNumOverride != -1) {
        stepNum = stepNumOverride;
    }

    if (stepNum == 1) {
        if (!startNewRecord || haventBeenHere) {
            doTrunc = true;
        }
	//  Add a time stamp
        ZZctml::addString(soln, "timeStamp", asctime(newtime));
	//  Add an identification XML element
	XML_Node* xmlID = eState_save_->writeIdentificationToXML();
	soln.mergeAsChild(*xmlID);
    }

    /*
     *  Add the globalTimeStep XML element with the global time step number as an attribute
     */
    ZZCantera::XML_Node& gts = soln.addChild("globalTimeStep");
    ii = ZZCantera::int2str(stepNum);
    gts.addAttribute("index", ii);
    std::string fmt = "%22.14E";

    /*
     *  Add the next time step's deltaT as a child XML element
     */
    ZZctml::addFloat(gts, "deltaTime_init_next", deltaTsubcycle_init_next_);
    /*
     *  Add this time step's initial deltaT try
     */
    ZZctml::addFloat(gts, "deltaTime_init_init", deltaTsubcycle_init_init_);
    //XML_Node& f = gts.addChild("deltaTime_init_next", deltaTsubcycle_init_next_ , fmt);
    //f.addAttribute("vtype", "float"); 
    /*
     *  Add the number of substep integrations as a child element
     */
    ZZctml::addInteger(gts, "numIntegrationSubCycles", numIntegrationSubCycles_final_final_);
    /*
     *  Add the child XML element, timeIncrement, to the globalTimeStep XML element.
     *  This includes the initial state, the final state, and the substep integration states.
     */
    if (!xmlTimeIncrementData_) {
        throw Electrode_Error("Electrode::writeSolutionTimeIncrement()",
                              "xmlTimeIncrementData_ is needed but wasn't created by the integrator");
    }
    gts.addChild(*xmlTimeIncrementData_);

    //       ------------------------- now finish up by writing the file -------------------------
    /*
     * Find the byte length of the current file
     */
    int length = 0;
    if (!doTrunc) {
        ifstream sr(fname.c_str());
        if (sr) {
            sr.seekg(0, ios::end);
            length = sr.tellg();
            sr.close();
        }
    }

    if (length > 0) {
        // Open the file and position it at the end of file
        fstream s(fname.c_str(), fstream::in | fstream::out | fstream::ate);

        if (startNewRecord) {
            /*
             * Backspace 8 characters
             */
            s.seekp(-8, ios_base::end);
            /*
             * Write a new electrodeOutput XML Element with an indentation of 2.
             */
            soln.write(s, 2);
        } else {
            //int pos = s.tellp();
            /*
             * Backspace 29 characters
             */
            s.seekp(-29, ios_base::end);
            /*
             * Write the globalTimeStep XML Element with an indentation of 4.
             */
            gts.write(s, 4);
            /*
             *  Now, write the end electrodeOutput
             */
            s.write("  </electrodeOutput>\n", 21);
        }
        /*
         *  Now, write the end ctml block to make the XML file whole again
         */
        s.write("</ctml>\n", 8);
        s.close();
    } else {
        /*
         * If the length is zero, write out the root XML element with all of the bells and whistles.
         */
        if (doTrunc) {
            fstream s(fname.c_str(), fstream::in | fstream::out | fstream::trunc);
            root.write(s);
            s.close();
        } else {
            fstream s(fname.c_str(), fstream::in | fstream::out | fstream::app);
            root.write(s);
            s.close();
        }
    }
    haventBeenHere = false;
}
//==================================================================================================================================
void Electrode::writeRestartFile(int stepNum, int solnNum, const std::string& nameAddition)
{

    //  solnNum is the global solution number. There can be more than one time dependent solution in an output file.

    if (!eState_save_) {
#ifdef DEBUG_MODE
        throw Electrode_Error("Electrode::writeRestartFile()", 
                              "Solution Output has been requested, "
                              "but it has never been set up by invoking electrode_stateSave_create()");
#endif
        return;
    }

    time_t aclock;
    ::time(&aclock); /* Get time in seconds */
    struct tm*   newtime = localtime(&aclock); /* Convert time to struct tm form */

    /*
     *  Determine the name of the solution file
     */
    std::string bname = "soln";
    if (baseNameSoln_ != "") {
        bname = baseNameSoln_ + "_Restart";
    }
    if (nameAddition != "") {
        bname += nameAddition;
    }

    /*
     *  Add in the domain number and cell number
     */
    std::string fname = (bname + "_" + int2str(electrodeDomainNumber_) + "_" + int2str(electrodeCellNumber_) + ".xml");
    
    //  Start a new root XML element and put the top elements in
    ZZCantera::XML_Node root("--");
    ZZCantera::XML_Node& ct = root.addChild("ctml");
    ZZCantera::XML_Node& soln = ct.addChild("electrodeOutput");

    /*
     *  Add the integer attribute, index, to the electrodeOutput element.
     *  -> so far we have an index of one. However, if there are continuation runs in the
     *     file and we want to add more than one
     */
    std::string ii = ZZCantera::int2str(solnNum);
    soln.addAttribute("index", ii);

    //  Add a time stamp
    ZZctml::addString(soln, "timeStamp", asctime(newtime));

   //   Add an identification XML element
#ifdef DEBUG_MODE
   XML_Node* xID = eState_save_->writeIdentificationToXML();
   soln.addChildToTree(xID);
#else
   soln.addChildToTree(eState_save_->writeIdentificationToXML());
#endif
    
    //  Add the globalTimeStep XML element with the global time step number as an attribute
    ZZCantera::XML_Node& gts = soln.addChild("globalTimeStep");
    ii = ZZCantera::int2str(stepNum);
    gts.addAttribute("index", ii);
    const std::string fmt = "%22.14E";

    //  Add this time step's first delta T attempt
    ZZctml::addFloat(gts, "deltaTime_init_init", deltaTsubcycle_init_init_);

    //  Add the next time step's deltaT as a child XML element
    ZZctml::addFloat(gts, "deltaTime_init_next", deltaTsubcycle_init_next_);

    //  Add the number of substep integrations as a child element
    ZZctml::addInteger(gts, "numIntegrationSubCycles", numIntegrationSubCycles_final_final_);
    /*
     *  Add the child XML element, timeIncrement, to the globalTimeStep XML element.
     *  This includes the initial state, the final state, and the substep integration states.
     */
    if (!xmlTimeIncrementData_) {
        throw Electrode_Error("Electrode::writeRestartFile()",
                              "xmlTimeIncrementData_ is needed but wasn't created by the integrator");
    }
    gts.addChild(*xmlTimeIncrementData_);

    //       ------------------------- now finish up by writing the file -------------------------

    std::fstream s(fname.c_str(), fstream::out | fstream::trunc);
    root.write(s);
    s.close();
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
