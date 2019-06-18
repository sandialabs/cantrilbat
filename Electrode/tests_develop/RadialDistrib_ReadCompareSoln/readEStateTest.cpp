/*
 * $Id: LiCoO2_Cathode_3_cc.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */
#include "gtest/gtest.h"

#include "zuzax/numerics/NonlinearSolver_JAC.h"
#include "zuzax/equilibrium.h"
#include "Electrode_Factory.h"
#include "importPL.h"
#include "BE_BlockEntry.h"
#include <iomanip>
#include "EState_XML.h"

using namespace std;
using namespace Zuzax;
using namespace esmodel;

// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;

void printUsage() {
    cout << "usage: readEstateTest [-h] [-help_cmdfile] [-d #] [mpequil.inp]"
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  electrodeCell.inp    : command file" << endl;
    cout << "                     : (if missing, assume mpequil.inp)" 
	 << endl;
    cout << endl;
}

//===================================================================================================
//! Unit test fixture for EState Comparison tests
class EStateComparisonTest : public testing::Test
{
public:
    //==============================================================================================================================
    EStateComparisonTest() :
	es1(0),
	es2(0),
        includeHist(true),
	printLvl(5),
	nDigits(6),
	molarAtol(1.0E-13),
	unitlessAtol(1.0E-13),
	time1(0.0),
        time2(0.0),
	eto1(0),
	eto2(0)
    {
    }
    //==============================================================================================================================
    ~EStateComparisonTest()
    {
	delete es1;
	delete es2;
	delete eto1;
	delete eto2;
    } 
    //==============================================================================================================================
    void shutdown()
    {
	delete es1;
	delete es2;
	delete eto1;
	delete eto2;
	es1 = es2 = 0;
	eto1 = eto2 = 0;
    } 
    //==============================================================================================================================
    bool ReadTwoFilesLastStepintoEState(const std::string& file1, const std::string& file2)
    {
	if (es1) {
	    delete es1;
	}
	es1 = readEState_XMLFile_LastStep("solnSaveA_0_0.xml", time1);
	if (!es1) {
	    return false;
	}
	if (es2) {
	    delete es2;
	}
	es2 = readEState_XMLFile_LastStep("solnSaveB_0_0.xml", time2);
	if (!es2) {
	    return false;
	}
	return true;
    }
    //==============================================================================================================================
    bool ReadFileIntoElectrodeTimeEvo(const std::string& file, int slot)
    {
	esmodel::ElectrodeTimeEvolutionOutput* eto = readFileElectrodeOutput(file);
	if (slot == 1) {
	    if (eto1) {
		delete eto1;
	    }
	    eto1 = eto;
	} else {
	    if (eto2) {
		delete eto2;
	    }
	    eto2 = eto;
	}
	return true;
    }
    //==============================================================================================================================
    bool WriteFilefromElectrodeTimeEvo(int slot, const std::string& file)
    {
	esmodel::ElectrodeTimeEvolutionOutput* eto = eto1;
	if (slot == 2) {
	    eto = eto2;
	}


	writeElectrodeOutputFile(file, *eto);

	return true;
    }
    //==============================================================================================================================
    bool CompareTwoEStates()
    {
	bool ok = es1->compareOtherState(es2, molarAtol, nDigits, includeHist, printLvl);
	return ok;
    }
    //==============================================================================================================================
    bool CompareTwoETimeEvolutions()
    {
	bool ok = eto1->compareOtherTimeEvolution(eto2, molarAtol, unitlessAtol, nDigits, includeHist, true, printLvl);
	return ok;
    }
    //==============================================================================================================================
    void setPrintLvl(int pp)
    {
	printLvl = pp;
    }
    //==============================================================================================================================
protected:


    Zuzax::EState* es1;
    Zuzax::EState* es2;
    bool includeHist;
    int printLvl;
    int nDigits;
    double molarAtol;
    double unitlessAtol;
    double time1;
    double time2;

    esmodel::ElectrodeTimeEvolutionOutput* eto1;
    esmodel::ElectrodeTimeEvolutionOutput* eto2;
};
//==================================================================================================================================
//==================================================================================================================================
TEST_F(EStateComparisonTest , ShouldWork)
{
    bool total_ok = true;
    setPrintLvl(0);
    bool r_ok = ReadTwoFilesLastStepintoEState("solnSaveA_0_0.xml", "solnSaveB_0_0.xml");
    bool ok = CompareTwoEStates();
  
    total_ok = r_ok && ok;
    EXPECT_TRUE(total_ok);

    shutdown();
}
//==================================================================================================================================
TEST_F(EStateComparisonTest , ReadWriteRead)
{
    bool total_ok = true;
    setPrintLvl(0);
    bool r_ok = ReadFileIntoElectrodeTimeEvo("solnSaveA_0_0.xml", 1);
    bool w_ok = WriteFilefromElectrodeTimeEvo(1, "solnSave_tmp.xml");
    bool r2_ok = ReadFileIntoElectrodeTimeEvo("solnSave_tmp.xml", 2);
    bool ok = CompareTwoETimeEvolutions();
    total_ok = r_ok && w_ok && r2_ok && ok;
    EXPECT_TRUE(total_ok);

    shutdown();
}
//==================================================================================================================================
TEST_F(EStateComparisonTest , Read2Fail)
{
    bool total_ok = true;
    setPrintLvl(0);
    bool r_ok = ReadFileIntoElectrodeTimeEvo("solnSaveA_0_0.xml", 1);

    bool r2_ok = ReadFileIntoElectrodeTimeEvo("solnSaveC_0_0.xml", 2);
    bool ok = CompareTwoETimeEvolutions();
    bool notOk = !ok;
    total_ok = r_ok && r2_ok && notOk;
    EXPECT_TRUE(total_ok);

    shutdown();
}
//==================================================================================================================================
int main(int argc, char **argv)
{

  //bool doCathode = false;
  string commandFileNet = "cell.inp";
  double time1, time2;
  string commandFileA = "anode.inp";
  // printed usage

  NonlinearSolver_JAC::s_TurnOffTiming = true;
  JacobianManager::s_print_Jac = true;


  Zuzax::XML_Node* xEout = getElectrodeOutputFile("solnSaveA_0_0.xml", 1);
  if (!xEout) {
      throw Electrode_Error("getElectrodeOutputFile", "Error");
  }
  std::string file1 = "solnSaveA_0_0.xml";
  std::string file2 = "solnSaveB_0_0.xml";

  Zuzax::EState* es1 = readEState_XMLFile_LastStep("solnSaveA_0_0.xml", time1);

  Zuzax::EState* es2 = readEState_XMLFile_LastStep("solnSaveB_0_0.xml", time2);
  bool includeHist = false;
  int printLvl = 5;
  int nDigits = 6;
  double molarAtol = 1.0E-13;
  printf(" comparison at times %g and %g\n", time1, time2);

  bool ok = es1->compareOtherState(es2, molarAtol, nDigits, includeHist, printLvl);
  delete es1;
  es1 = 0;
  delete es2;
  es2 = 0;


  if (!ok) {
      printf("WARNING files are different\n");
  }

  EState_ID_struct e_id;
  e_id.readIdentificationFromXML(*xEout);

  EState* es = newEStateObject(e_id.EState_Type_String);
  es->readIdentificationFromXML(*xEout); 

  int globalTimeStepNum = 0;
  Zuzax::XML_Node* x = selectLastGlobalTimeStepInterval(xEout, globalTimeStepNum);

  double timeVal;
  Zuzax::XML_Node* xSt = locateTimeLast_GlobalTimeStepIntervalFromXML(*x, timeVal, 1);
  es->readStateFromXML(*xSt);

  delete es; 
  es = 0;
  printf("read global time step num %d\n", globalTimeStepNum);

  //
  // Enter google test environment
  //
  testing::InitGoogleTest(&argc, argv);
  int res = RUN_ALL_TESTS();
  
  Zuzax::appdelete();


  return res;

} 
