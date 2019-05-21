/**
 *  @file changeHf.cpp
 *
 */

//  Example 
//
//  Read an XML file. 
//  Look for the reference state thermodynamics of a single
//  species and then read it in.
//  Then, change the 298.15 Hf value
//  Then, print out the xml file with the changed 298.15 HF value
//  reporting what has happened.
//


#include "LE_PickList.h"
#include "BlockEntry.h"


#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "changeHf.h"

#include "zuzax/kernel/GeneralSpeciesThermo.h"
#include "zuzax/kernel/SpeciesThermoFactory.h"
#include "zuzax/kernel/Mu0Poly.h"
#include "zuzax/base/ctml.h"
#include "zuzax/thermo/SpeciesThermo.h"

using namespace Cantera;
using namespace std;
using namespace ctml;
using namespace BEInput;

int DebugPrinting = 0;
#ifdef DEBUG_HKM
int iDebug_HKM = 0;
#endif


UnitsIO UIO;

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
static void printUsage()
{
    printf("changeHf Usage:\n");
    printf("\t changeHf [file]\n");
    printf("\t\twhere the default file name is changeHf.inp\n");
    printf("\n");
}


/** 
 * Install a NASA polynomial thermodynamic property
 * parameterization for species k into a SpeciesThermo instance.
 *
 * There may be 1 or 2 regions. If f1ptr == 0 there is only 
 * one region.
 */
static void installNasa1ThermoFromXML(std::string speciesName,
				      SpeciesThermo& sp, int k, 
				      const XML_Node* f0ptr) {
    doublereal tmin, tmax;
    const XML_Node& f0 = *f0ptr;
    tmin = fpValue(f0["Tmin"]);
    tmax = fpValue(f0["Tmax"]);
 
    vector_fp c0;
    ctml::getFloatArray(f0.child("floatArray"), c0, false);
    vector_fp c(7);
    doublereal p0 = OneAtm;
    c[0] = c0[5];
    c[1] = c0[6];
    copy(c0.begin(), c0.begin()+5, c.begin() + 2);
    sp.install(speciesName, k, NASA1, &c[0], tmin, tmax, p0);
}

static void rewriteNasa1ThermoToXML(XML_Node* f0ptr, double *c) {
    doublereal tmin, tmax;
    XML_Node& f0 = *f0ptr;
    XML_Node *nasaXMLptr = f0.findByName("NASA");
    XML_Node& nasaXML = *nasaXMLptr;
    tmin = fpValue(nasaXML["Tmin"]);
    tmax = fpValue(nasaXML["Tmax"]);
    XML_Node &faXML = nasaXML.child("floatArray"); 
    nasaXML.removeChild(&faXML);

    vector_fp c0(7);
    c0[5] = c[0];
    c0[6] = c[1];
    c0[0] = c[2];
    c0[1] = c[3];
    c0[2] = c[4];
    c0[3] = c[5];
    c0[4] = c[6];

    string title = "coeffs";
    addFloatArray(nasaXML, title, 7, &c0[0]);
 
}


 /** 
  * Install a NASA polynomial thermodynamic property
  * parameterization for species k into a SpeciesThermo instance.
  *
  * There may be 1 or 2 regions. If f1ptr == 0 there is only 
  * one region.
  */
static void installNasa2ThermoFromXML(string speciesName,
				      SpeciesThermo& sp, int k, 
				      const XML_Node* f0ptr, 
				      const XML_Node* f1ptr) {
    doublereal tmin0, tmax0, tmin1, tmax1, tmin, tmid, tmax;

    const XML_Node& f0 = *f0ptr;
    if (!f1ptr) {
      printf("DuelRange requrired\n");
      exit(-1);
    }
    bool dualRange = false;
    if (f1ptr) {dualRange = true;}
    tmin0 = fpValue(f0["Tmin"]);
    tmax0 = fpValue(f0["Tmax"]);
    tmin1 = tmax0;
    tmax1 = tmin1 + 0.0001;
    if (dualRange) {
      tmin1 = fpValue((*f1ptr)["Tmin"]);
      tmax1 = fpValue((*f1ptr)["Tmax"]);
    }

    vector_fp c0, c1;
    if (fabs(tmax0 - tmin1) < 0.01) {
      tmin = tmin0;
      tmid = tmax0;
      tmax = tmax1;
      getFloatArray(f0.child("floatArray"), c0, false);
      if (dualRange)
	  getFloatArray(f1ptr->child("floatArray"), c1, false);
      else {
	c1.resize(7,0.0);
	copy(c0.begin(), c0.end(), c1.begin());
      }
    }
    else if (fabs(tmax1 - tmin0) < 0.01) {
      tmin = tmin1;
      tmid = tmax1;
      tmax = tmax0;
      getFloatArray(f1ptr->child("floatArray"), c0, false);
      getFloatArray(f0.child("floatArray"), c1, false);
    }
    else {
      throw ZuzaxError("installNasaThermo",
			 "non-continuous temperature ranges.");
    }
    vector_fp c(15);
    c[0] = tmid;
    doublereal p0 = OneAtm;
    c[1] = c0[5];
    c[2] = c0[6];
    copy(c0.begin(), c0.begin()+5, c.begin() + 3);
    c[8] = c1[5];
    c[9] = c1[6];
    copy(c1.begin(), c1.begin()+5, c.begin() + 10);
    sp.install(speciesName, k, NASA2, &c[0], tmin, tmax, p0);
}

#ifdef INCL_NASA96

/** 
 * Install a NASA96 polynomial thermodynamic property
 * parameterization for species k into a SpeciesThermo instance.
 */
static void installNasa96ThermoFromXML(string speciesName,
				       SpeciesThermo& sp, int k, 
				       const XML_Node* f0ptr, 
				       const XML_Node* f1ptr) {
    doublereal tmin0, tmax0, tmin1, tmax1, tmin, tmid, tmax;

    const XML_Node& f0 = *f0ptr;
    bool dualRange = false;
    if (f1ptr) {dualRange = true;}
    tmin0 = fpValue(f0["Tmin"]);
    tmax0 = fpValue(f0["Tmax"]);
    tmin1 = tmax0;
    tmax1 = tmin1 + 0.0001;
    if (dualRange) {
      tmin1 = fpValue((*f1ptr)["Tmin"]);
      tmax1 = fpValue((*f1ptr)["Tmax"]);
    }

    vector_fp c0, c1;
    if (fabs(tmax0 - tmin1) < 0.01) {
      tmin = tmin0;
      tmid = tmax0;
      tmax = tmax1;
      getFloatArray(f0.child("floatArray"), c0, false);
      if (dualRange)
	  getFloatArray(f1ptr->child("floatArray"), c1, false);
      else {
	c1.resize(7,0.0);
	copy(c0.begin(), c0.end(), c1.begin());
      }
    }
    else if (fabs(tmax1 - tmin0) < 0.01) {
      tmin = tmin1;
      tmid = tmax1;
      tmax = tmax0;
      getFloatArray(f1ptr->child("floatArray"), c0, false);
      getFloatArray(f0.child("floatArray"), c1, false);
    }
    else {
      throw ZuzaxError("installNasaThermo",
			 "non-continuous temperature ranges.");
    }
    array_fp c(15);
    c[0] = tmid;
    doublereal p0 = OneAtm;
    c[1] = c0[5];
    c[2] = c0[6];
    copy(c0.begin(), c0.begin()+5, c.begin() + 3);
    c[8] = c1[5];
    c[9] = c1[6];
    copy(c1.begin(), c1.begin()+5, c.begin() + 10);
    sp.install(speciesName, k, NASA, c.begin(), tmin, tmax, p0);
}

#endif


/** 
 * Install a Shomate polynomial thermodynamic property
 * parameterization for species k.
 */
static void installShomateThermoFromXML(string speciesName, 
					SpeciesThermo& sp, int k, 
					const XML_Node* f0ptr,
					const XML_Node* f1ptr) {
    doublereal tmin0, tmax0, tmin1, tmax1, tmin, tmid, tmax;

    const XML_Node& f0 = *f0ptr;
    bool dualRange = false;
    if (f1ptr) {dualRange = true;}
    tmin0 = fpValue(f0["Tmin"]);
    tmax0 = fpValue(f0["Tmax"]);
    tmin1 = tmax0;
    tmax1 = tmin1 + 0.0001;
    if (dualRange) {
      tmin1 = fpValue((*f1ptr)["Tmin"]);
      tmax1 = fpValue((*f1ptr)["Tmax"]);
    }

    vector_fp c0, c1;
    if (fabs(tmax0 - tmin1) < 0.01) {
      tmin = tmin0;
      tmid = tmax0;
      tmax = tmax1;
      getFloatArray(f0.child("floatArray"), c0, false);
      if (dualRange)
	  getFloatArray(f1ptr->child("floatArray"), c1, false);
      else {
	c1.resize(7,0.0);
	copy(c0.begin(), c0.begin()+7, c1.begin());
      }
    }
    else if (fabs(tmax1 - tmin0) < 0.01) {
      tmin = tmin1;
      tmid = tmax1;
      tmax = tmax0;
      getFloatArray(f1ptr->child("floatArray"), c0, false);
      getFloatArray(f0.child("floatArray"), c1, false);
    }
    else {
      throw ZuzaxError("installShomateThermo",
			 "non-continuous temperature ranges.");
    }
    vector_fp c(15);
    c[0] = tmid;
    doublereal p0 = OneAtm;
    copy(c0.begin(), c0.begin()+7, c.begin() + 1);
    copy(c1.begin(), c1.begin()+7, c.begin() + 8);
    sp.install(speciesName, k, SHOMATE, &c[0], tmin, tmax, p0);
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/


/** 
 * Install a constant-cp thermodynamic property
 * parameterization for species k.
 */
static void installSimpleThermoFromXML(string speciesName, 
				       SpeciesThermo& sp, int k, 
				       const XML_Node& f) {
    doublereal tmin, tmax;
    tmin = fpValue(f["Tmin"]);
    tmax = fpValue(f["Tmax"]);
    if (tmax == 0.0) tmax = 1.0e30;

    vector_fp c(4);
    c[0] = getFloat(f, "t0", "-");
    c[1] = getFloat(f, "h0", "-");
    c[2] = getFloat(f, "s0", "-");
    c[3] = getFloat(f, "cp0", "-");
    doublereal p0 = OneAtm;
    sp.install(speciesName, k, SIMPLE, &c[0], tmin, tmax, p0);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
 *
 */
void installThermoForSpecies(int k, const XML_Node& s, 
			     SpeciesThermo& spthermo) {
    /*
     * Check to see that the species block has a thermo block
     * before processing. Throw an error if not there.
     */
    if (!(s.hasChild("thermo"))) {
      throw UnknownSpeciesThermoModel("installSpecies", 
				      s["name"], "<nonexistent>");
    }
    const XML_Node& thermo = s.child("thermo");
    const vector<XML_Node*>& tp = thermo.children();
    int nc = static_cast<int>(tp.size());
    if (nc == 1) {
      const XML_Node* f = tp[0];
      if (f->name() == "Shomate") {
	installShomateThermoFromXML(s["name"], spthermo, k, f, 0);
      }
      else if (f->name() == "const_cp") {
	installSimpleThermoFromXML(s["name"], spthermo, k, *f);
      }
      else if (f->name() == "NASA") {
	installNasa1ThermoFromXML(s["name"], spthermo, k, f);
      }
      else if (f->name() == "Mu0") {
	installMu0ThermoFromXML(s["name"], spthermo, k, f);
      }
      else {
	throw UnknownSpeciesThermoModel("installSpecies", 
					s["name"], f->name());
      }
    }
    else if (nc == 2) {
      const XML_Node* f0 = tp[0];
      const XML_Node* f1 = tp[1];
      if (f0->name() == "NASA" && f1->name() == "NASA") {
	installNasa2ThermoFromXML(s["name"], spthermo, k, f0, f1);
      } 
      else if (f0->name() == "Shomate" && f1->name() == "Shomate") {
	installShomateThermoFromXML(s["name"], spthermo, k, f0, f1);
      } 
      else {
	throw UnknownSpeciesThermoModel("installSpecies", s["name"], 
					f0->name() + " and "
					+ f1->name());
      }
    }
    else {
      throw UnknownSpeciesThermoModel("installSpecies", s["name"], 
				      "multiple");
    }
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
 *
 */
void changeH_RT(GeneralSpeciesThermo *gsp, XML_Node* s, 
		double OldH_RT, double NewH_RT) {

    vector_fp c(50);
    int type;
    double minTemp, maxTemp, refPressure;
    gsp->reportParams(0, type, &c[0],
		      minTemp, maxTemp, refPressure);

    if (type == NASA1) {
      printf("NASA1\n");
      c[0] = c[0] + (NewH_RT - OldH_RT) * 298.15;
      gsp->install(IOO.SpeciesName, 1, NASA1, &c[0], minTemp, maxTemp,
		   refPressure);  
    }

    double RT = GasConstant * 298.15;
    double H_RT[2], Cp_R[2], S_R[2];
    gsp->update(298.15, (double *) &Cp_R, (double *)&H_RT, (double *)&S_R);

    printf("Hold = %13.8g J/kmol\n", H_RT[0] * RT);
    printf("Hnew = %13.8g J/kmol\n", H_RT[1] * RT);

    gsp->reportParams(1, type, &c[0],
		      minTemp, maxTemp, refPressure);
    if (type == NASA1) {
      rewriteNasa1ThermoToXML( s, &c[0]);
    } else {
      printf("Case not handled yet\n");
      exit(-1);
    }
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

int main(int argc, char** argv) {
    std::string infile = "changeHf.inp";  
#ifdef DEBUG_HKM
    DebugPrinting = true;
#endif
    // look for command-line options
    if (argc > 1) {
      string tok;
      for (int j = 1; j < argc; j++) {
	tok = string(argv[j]);
	if (tok[0] == '-') {
	  int nopt = tok.size();
	  for (int n = 1; n < nopt; n++) {
	    if (tok[n] == 'h') {
	      printUsage();
	      exit(0);
	    } else {
	      printUsage();
	      exit(1);
	    }
	  }
	} else if (infile == "" || infile == "changeHf.inp") {
	  infile = tok;
	}
	else {
	  printUsage();
	  exit(1);
	}
      }
    }
 
    /*
     * Set up and process the file
     */
    FILE *fff = fopen(infile.c_str(), "r");
    BlockEntry *cf = new BlockEntry("command_file");
    setup_input(cf);
    process_input(cf, fff);
    delete cf;
    fclose (fff);


    XML_Node *xc = new XML_Node();
    const char *fn = IOO.FileName.c_str();
    string path = findInputFile(fn);
    ifstream fin(path.c_str());
    if (!fin) {
      throw ZuzaxError("changeHf","could not open "
			 +path+" for reading.");
    }
    /*
     * Make a complete copy of the xml file
     */
    xc->build(fin);
    fin.close();
    XML_Node *xd = new XML_Node();
    xc->copy(xd);

    /*
     * Delete the original copy
     */
    delete xc;
    xc = 0;

    // Find the species in the database by name.
    XML_Node* s = xd->findByAttr("name", IOO.SpeciesName);
    if (s) {
      if (DebugPrinting) {
	printf("Found species %s\n", IOO.SpeciesName.c_str());
      }
    } else {
      printf("No species %s found\n", IOO.SpeciesName.c_str()); 
      exit(-1);
    }
    const XML_Node *sc = s;
    GeneralSpeciesThermo *gsp = new GeneralSpeciesThermo();
  

    installThermoForSpecies(0, *sc, (SpeciesThermo &) *gsp);


    double RT = GasConstant * 298.15;
    double H_RT, Cp_R, S_R;
    gsp->update(298.15, &Cp_R, &H_RT, &S_R);

    /*
     * Print out the old value of the Heat of Formation
     */
    if (DebugPrinting) {
      printf("Found Hfold = %13.8g J/kmol\n", H_RT * RT);
    }
  
    /*
     * Calculate the New value of the Heat of Formation
     */
    double NewH_RT = H_RT;
    if (IOO.DeltaValue != 0.0) {
      NewH_RT = (IOO.DeltaValue + H_RT * RT) / RT;
    } else {
      NewH_RT = (IOO.HfValue) / RT;
    }

    /*
     * Print out the new value of the Heat of Formation
     */
    if (DebugPrinting) {
      printf("Setting Hfnew = %13.8g J/kmol\n", NewH_RT * RT);
    }

    /*
     * Change the Heat of Formation
     */
    changeH_RT(gsp, s, H_RT, NewH_RT);

    /*
     * Write the revised XML file out to the destination file.
     */
    ofstream tout;
    const char *destfilename = IOO.DestFileName.c_str();
    tout.open(destfilename);
    xd->write(tout);
    tout.close();

    /*
     * Clean up
     */
    delete xd;
    appdelete();

    return 0;
}
/***********************************************************/
