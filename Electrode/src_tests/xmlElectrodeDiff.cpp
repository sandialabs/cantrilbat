/**
 * @file xmlElectrodeDiff.cpp
 *  Main routine for carrying out an xml diff on saved electrode files in esmodel Electrode XML format.
 *  (see \ref electrode_xml_format ).
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

/*
 *
 *  xmlSolnDiff File1.xml File2.xml
 *
 *  Compares the variable values in two Cantera solution xml 
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


#include "Electrode.h"


#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include <fstream>
#include <unistd.h>

#include "EState.h"
#include "EState_XML.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
using namespace esmodel;


#if defined(__CYGWIN__)
#include <getopt.h>
#endif


#include "tok_input_util.h"

#ifdef WIN32
#pragma warning(disable:4996)
#endif
#ifndef MAX
#  define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#  define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

int Debug_Flag = 1;
double grtol = 1.0E-6;
double gatol = 1.0E-12;
std::map<std::string, double> VarRtol;
std::map<std::string, double> VarAtol;

//======================================================================================================================
#ifdef WINMSVC
/*
 * Windows doesn't have getopt(). This is an incomplete version that
 * does enough to handle required functionality.
 */
int optind = -1;
char *optarg = 0;

int getopt(int argc, char **argv, const char *) {
    static int currArg = 1;
    static int currOptInd = 1;
    string tok;
    static int charPos = 0;
    int rc = -1;
    if (currArg >= argc) {
      optarg = 0;
      return -rc;
    }
    tok = string(argv[currArg]);
    currOptInd = currArg+1;
    if (currOptInd > argc - 1) {
      currOptInd = -1;
      optarg = 0;
    } else {
      optarg = argv[currArg+1];
    }
    size_t len = strlen(tok.c_str());
    if (charPos == 0) {
      bool found = false;
      do {
	tok = string(argv[currArg]);
	len = strlen(tok.c_str());
	if (len > 1 && tok[0] == '-') {
	  found = true;
	  charPos = 1;
	  if (len > 2 && tok[1] == '-') {
	    charPos = 2;
	  }
	} else {
	  if (optind == -1) {
	    optind = currArg;
	  }
	}
	if (!found) {
	  if (currArg < (argc-1)) {
	    currArg++;
	  } else {
            optarg = 0;
	    return -1;
	  }
	}
      } while (!found);
    }
 
    rc = tok[charPos];
    if (charPos < static_cast<int>(len - 1)) {
      charPos++;
    } else {
      charPos = 0;
    }
    return rc;
}

#endif


//======================================================================================================================
XML_Node *readXML(std::string inputFile) {

  if (inputFile.size() == 0) {
    throw CanteraError("constructXML",  "input file is null");
  }
  string path = findInputFile(inputFile);
  std::ifstream fin(path.c_str());
  if (!fin) {
    throw CanteraError("HMWSoln:constructPhaseFile","could not open "
                       +path+" for reading.");
  } 
  /*
   * The phase object automatically constructs an XML object.
   * Use this object to store information.
   */
  XML_Node *fxml = new XML_Node();
  fxml->build(fin);
  return fxml;
}

//======================================================================================================================
static void print_usage() {
  printf("\t\n");
  printf("  xmlElectrodeDiff [-h] [-a atol] [-r rtol] File1.xml File2.xml\n");
  printf("\t\n");
  printf("\tCompares two solution files in esmodel Electrode XML format.\n");
  printf("\tThe comparison is done using a weighted norm basis.\n");
  printf("\t\n");
  printf("\tThe two files should be basically equal. However, File1.xml is\n");
  printf("\ttaken as the reference file that has precedence, when there is\n");
  printf("\tsomething to be decided upon.\n");
  printf("\t\n");
  printf("\t Arguments:\n");
  printf("\t  -h      = Usage info\n");
  printf("\t  -d      = printLvl \n");
  printf("\t  -a atol = Set absolute tolerance parameter - default = 1.0E-12\n");
  printf("\t  -r rtol = Set relative tolerance parameter - default = 1.0E-6\n");
  printf("\t\n");
  printf("\t Shell Return Values:\n");
  printf("\t   0 = Comparison was successful\n");
  printf("\t   1 = One or more time intervals failed the comparison\n");
  printf("\t\n");
}
//====================================================================================================================
int main(int argc, char *argv[])
{
  int opt_let;
  char  *fileName1=NULL, *fileName2=NULL;
  XML_Node *Xfp1 = 0;
  int printLvl = 0;
  XML_Node *Xfp2 = 0;
  double atol_arg = 0.0, rtol_arg = 0.0;
  int id = 0;
  int id2 = 0;
  char *ggg = 0;
  char *rrr = 0;
  /*
   * Interpret command line arguments
   */
  /* Loop over each command line option */
  while((opt_let = getopt(argc, argv, "ha:r:d:")) != EOF) {

    /* case over the option letter */
    switch(opt_let) {
      
    case 'h':
      /* Usage info was requested */
      print_usage();
      exit(0);

    case 'a':
      /* atol parameter */
 
      ggg = optarg;
      //printf("a = %s\n", ggg);
      id = sscanf(ggg,"%lg", &atol_arg);
      if (id != 1) {
	printf(" atol param bad: %s\n", ggg);
	exit(-1);
      }
      gatol = atol_arg;
      break;

    case 'r':
      /* rtol parameter */
 
      rrr = optarg;
      //printf("r = %s\n", ggg);
      id2 = sscanf(rrr,"%lg", &rtol_arg);
      if (id2 != 1) {
	printf(" rtol param bad: %s\n", rrr);
	exit(-1);
      }
      grtol = rtol_arg;
      break;


    case 'd':
      /* printlvl parameter */
 
      rrr = optarg;
      id2 = sscanf(rrr,"%lg", &rtol_arg);
      if (id2 != 1) {
	printf(" rtol param bad: %s\n", rrr);
	exit(-1);
      }
      printLvl = rtol_arg;
      break;

    
      
    default:
      /* Default case. Error on unknown argument. */
      printf("default called opt_let = %c\n", opt_let);
      fprintf(stderr, "ERROR in command line usuage:\n");
      print_usage();
      return 0; 
    } /* End "switch(opt_let)" */

  } /* End "while((opt_let=getopt(argc, argv, "i")) != EOF)" */

  if (optind !=  argc-2) {
    print_usage();
    exit(-1);
  } else {
    fileName1 = argv[argc-2];
    fileName2 = argv[argc-1];
  }

  /*
   *      Print Out Header
   */
  printf("\n");
  printf("----------------------------------------------------------------------\n");
  printf("xmlElectrodeDiff: XML Soln File comparison utility program\n");
  printf("             Version \n");
  printf("             Harry K. Moffat Div. 1516 Sandia National Labs\n");
  printf("           \n");
  printf("             First  Cantera XML Electrode File = %s\n", fileName1);
  printf("             Second Cantera XML Electrode File = %s\n", fileName2); 
  printf("\n");
  printf("             Absolute tol = %g\n", gatol);
  printf("             Relative tol = %g\n", grtol);
  printf("----------------------------------------------------------------------\n");
  printf("\n");  
  
  /*
   *  Open up the two xml Files #1 and #2
   */

  if (!(Xfp1 = readXML(fileName1))) {
    fprintf(stderr,"Error opening up file1, %s\n", fileName1);
    exit(-1);
  }
  if (!(Xfp2 = readXML(fileName2))) {
    fprintf(stderr, "Error opening up file2, %s\n", fileName2);
    exit(-1);
  }
  esmodel::ElectrodeTimeEvolutionOutput* eto1 = 0;
  esmodel::ElectrodeTimeEvolutionOutput* eto2 = 0; 

  try {
      eto1 = readXMLElectrodeOutput(*Xfp1);
  } catch (CanteraError) {
    showErrors();
    printf("xmlElectrodeDiff ERROR:: Error while reading first file, %s. Exiting program\n", fileName1);
    return -1;
  }
  try {
      eto2 = readXMLElectrodeOutput(*Xfp2);
  } catch (CanteraError) {
    showErrors();
    printf("xmlElectrodeDiff ERROR:: Error while reading second file, %s. Exiting program\n", fileName2);
    return -1;
  }

  if (eto1->numGlobalTimeIntervals_ != eto2->numGlobalTimeIntervals_) {
      printf("Warning: number of time intervals are unequal\n");
  }
  double molarAtol = 1.0E-20;
  double digits = -log10(grtol);
  int nDigits = digits;
  int includeHist = 1;
  int compareType = 1;
 
  int numZonesNeededToPass = std::min(eto1->numGlobalTimeIntervals_, eto2->numGlobalTimeIntervals_);
  int numPassed =  numZonesNeededToPass;
  bool ok;
  try {
      ok = eto1->compareOtherTimeEvolutionSub(eto2, numPassed, molarAtol, gatol, nDigits,includeHist, compareType, printLvl);
  } catch (CanteraError) {
    showErrors();
    printf("xmlElectrodeDiff ERROR:: Error while comparing first solution %s to second solution %s. Exiting program\n",
	   fileName1, fileName2);
    return -1;
  }

  if (ok) {
      printf("xmlElectrodeDiff: Passed, a total of %d zones were the same\n",  numZonesNeededToPass);
  } else {
      int numDiff = numZonesNeededToPass - numPassed;
      printf("xmlElectrodeDiff: Failed, a total of %d zones were the same, a total of %d intervals were different\n",
	     numZonesNeededToPass, numDiff);
  }
  int iReturn = !ok;
  if (ok == true) {
      iReturn = 0;
  } else {
      iReturn = 1;
  }
 
  return iReturn;

}
//======================================================================================================================
