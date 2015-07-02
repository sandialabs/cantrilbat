/*====================================================================
 * ------------------------
 * | CVS File Information | 
 * ------------------------
 * $RCSfile: xmlSolnDiff.cpp,v $
 * $Author: hkmoffa $
 * $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 * $Revision: 5 $
 * $Name:  $
 *====================================================================*/
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


using namespace std;
using namespace Cantera;



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
double grtol = 1.0E-3;
double gatol = 1.0E-9;
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
double rtolVar(std::string var) {
  std::map<std::string, double>::const_iterator i = VarRtol.find(var);
  if (i != VarRtol.end()) {
    return i->second;
  } else {
    VarRtol[var] = grtol;
  }
  return  VarRtol[var];
}
//======================================================================================================================
double atolVar(std::string var) {
  std::map<std::string, double>::const_iterator i = VarAtol.find(var);
  if (i != VarAtol.end()) {
    return i->second;
  } else {
    VarAtol[var] = grtol;
  }
  return  VarAtol[var];
}
//======================================================================================================================
int diff_double(double d1, double d2, double rtol, double atol)

/*
 * Compares 2 doubles. If they are not within tolerance, then this
 * function returns true. 
 */
{
  if (fabs(d1-d2) > (atol + rtol * 0.5 * (fabs(d1) + fabs(d2)))) return 1;
  return 0;
}
//======================================================================================================================
static int diff_double_slope(double d1, double d2, double rtol, 
			     double atol, double xtol, double slope1, double slope2)

/*
 * Compares 2 doubles. If they are not within tolerance, then this
 * function returns true. 
 */
{
  double atol2 = xtol*(fabs(slope1) + fabs(slope2));
  if (fabs(d1-d2) > (atol + atol2 + rtol * 0.5 * (fabs(d1) + fabs(d2)))) return 1;
  return 0;
}
//======================================================================================================================
static double calc_rdiff(double d1, double d2, double rtol, double atol)
{
  double rhs, lhs;
  rhs = fabs(d1-d2);
  lhs = atol + rtol * 0.5 * (fabs(d1) + fabs(d2));
  return (rhs/lhs);
}
//======================================================================================================================
static double get_atol(std::vector<double> &values, const int nvals,const double atol) {
  double sum = 0.0, retn;
  if (nvals <= 0) return gatol;
  for (int i = 0; i < nvals; i++) {
    retn = values[i];
    sum += retn * retn;
  }
  sum /= nvals;
  retn = sqrt(sum);
  return ((retn + 1.0) * atol);
}
//======================================================================================================================
bool compareFieldVectorFlts(std::vector<double> & comp1,
			    std::vector<double> & comp2,
			    std::vector<double> & xpos,
			    std::string varName) {

  double rtol = rtolVar(varName);
  double atol = atolVar(varName);
 
  double max_diff = 0.0;
  int ndiff = 0;
  double rel_diff = 0.0;
  double xatol, slope1, slope2;

  int nDataRows1 = comp1.size();
  int nDataRows2 = comp2.size();
  int nDataRowsMIN = MIN(nDataRows1, nDataRows2);
  int nDataRowsMAX = MAX(nDataRows1, nDataRows2);
  // Adjust for large values
  double atol_j =  get_atol(comp1, nDataRows1, atol);
  atol_j = MIN(atol_j, get_atol(comp2, nDataRows2, atol));
  int jmax = 0;

  for (int j = 0; j < nDataRowsMIN; j++) {
    
    slope1 = 0.0;
    slope2 = 0.0;

    if (j == 0) {
      slope1 = (comp1[j+1] - comp1[j])/(xpos[j+1] - xpos[j]);
      slope2 = (comp2[j+1] - comp2[j])/(xpos[j+1] - xpos[j]);
      xatol = fabs(grtol *(xpos[1] - xpos[0]));
    } else if (j == (nDataRowsMIN-1)) {
      slope1 = (comp1[j] - comp1[j-1])/(xpos[j] - xpos[j-1]);
      slope2 = (comp2[j] - comp2[j-1])/(xpos[j] - xpos[j-1]);
      xatol = fabs(grtol *(xpos[j] - xpos[j-1]));
    } else {
      slope1 = (comp1[j+1] - comp1[j-1])/(xpos[j+1] - xpos[j-1]);
      slope2 = (comp2[j+1] - comp2[j-1])/(xpos[j+1] - xpos[j-1]);
      xatol = fabs(grtol * 0.5 * (xpos[j+1] - xpos[j-1]));
    }
    bool notOK = diff_double_slope(comp1[j], comp2[j], 
			      rtol, atol_j, xatol, slope1, slope2);
    if (notOK) {
      ndiff++;
      rel_diff = calc_rdiff((double) comp1[j], (double) comp2[j], rtol, atol_j);
      if (rel_diff > max_diff) {
	jmax = j;
	max_diff = rel_diff;
      }
      if (ndiff < 10 || (jmax == j)) {
	printf("\tfield variable %s at Node  %d ", varName.c_str(), j);
	printf(" differs: %g %g\n", comp1[j],  comp2[j]);
      }
    }
  }
    
  if (nDataRowsMIN != nDataRowsMAX) {
    ndiff += nDataRowsMAX - nDataRowsMIN;
    if (ndiff < 10) {
      if (nDataRows1 > nDataRows2) {
	for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
	  printf("\tColumn variable %s at data row %d ", varName.c_str(), j + 1);
	  printf(" differ: %g      NA\n", comp1[j]);
	}
      } else {
	for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
	  printf("\tColumn variable %s at data row %d ", varName.c_str(), j + 1);
	  printf(" differ: NA     %g \n", comp2[j]);
	}
      }
    }
  }
  if (ndiff == 0) {
    printf("Field Variable %s passed test on %d Nodes\n",  varName.c_str(), nDataRowsMIN);
  } else {
    printf("Field Variable %s failed test on %d Nodes out of %d Nodes\n",  varName.c_str(), ndiff, nDataRowsMIN);
  }
  return (ndiff == 0);
}
//======================================================================================================================
bool compareVectorFlts(std::vector<double> & comp1,
		       std::vector<double> & comp2,
		       std::string varName) {

  double rtol = rtolVar(varName);
  double atol = atolVar(varName);
 
  double max_diff = 0.0;
  int ndiff = 0;
  double rel_diff = 0.0;

  int nDataRows1 = comp1.size();
  int nDataRows2 = comp2.size();
  int nDataRowsMIN = MIN(nDataRows1, nDataRows2);
  int nDataRowsMAX = MAX(nDataRows1, nDataRows2);
  // Adjust for large values
  double atol_j =  get_atol(comp1, nDataRows1, atol);
  atol_j = MIN(atol_j, get_atol(comp2, nDataRows2, atol));
  int jmax = 0;

  for (int j = 0; j < nDataRowsMIN; j++) {
    
 
    bool notOK = diff_double(comp1[j], comp2[j], rtol, atol_j);
    if (notOK) {
      ndiff++;
      rel_diff = calc_rdiff((double) comp1[j], (double) comp2[j], rtol, atol_j);
      if (rel_diff > max_diff) {
	jmax = j;
	max_diff = rel_diff;
      }
      if (ndiff < 10 || (jmax == j)) {
	printf("\tVariable %s at Node  %d ", varName.c_str(), j);
	printf(" differs: %g %g\n", comp1[j],  comp2[j]);
      }
    }
  }
    
  if (nDataRowsMIN != nDataRowsMAX) {
    ndiff += nDataRowsMAX - nDataRowsMIN;
    if (ndiff < 10) {
      if (nDataRows1 > nDataRows2) {
	for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
	  printf("\tvariable %s at data row %d ", varName.c_str(), j + 1);
	  printf(" differ: %g      NA\n", comp1[j]);
	}
      } else {
	for (int j = nDataRowsMIN; j < nDataRowsMAX; j++) {
	  printf("\tvariable %s at data row %d ", varName.c_str(), j + 1);
	  printf(" differ: NA     %g \n", comp2[j]);
	}
      }
    }
  }
  if (ndiff == 0) {
    printf("Variable %s passed test on %d Nodes\n", varName.c_str(), nDataRowsMIN);
  } else {
    printf("Variable %s failed test on %d Nodes out of %d Nodes\n", varName.c_str(), ndiff, nDataRowsMIN);
  }
  return (ndiff == 0);
}

//======================================================================================================================
XML_Node *getSimul(XML_Node *xmlTop, std::string id_tag) {
  XML_Node *xctml = xmlTop;
  if (xctml->name() != "ctml") {
    xctml = xmlTop->findByName("ctml");
    if (!xctml) {
      throw CanteraError("countSimulations","can't find ctml node");
    }
  }
  XML_Node* node = xctml->findNameID("simulation", id_tag);
  return node;
}

//======================================================================================================================
int countSimulations(XML_Node *xmlTop) {
  XML_Node *xctml = xmlTop;
  if (xctml->name() != "ctml") {
    xctml = xmlTop->findByName("ctml");
    if (!xctml) {
      throw CanteraError("countSimulations","can't find ctml node");
    }
  }
  std::vector<XML_Node*> ccc;
  xctml->getChildren("simulation", ccc);
  int sz = ccc.size();
  return sz;
}   
//======================================================================================================================
XML_Node *getSimulNum(XML_Node *xmlTop, int num) {
  XML_Node *xctml = xmlTop;
  if (xctml->name() != "ctml") {
    xctml = xmlTop->findByName("ctml");
    if (!xctml) {
      throw CanteraError("countSimulations","can't find ctml node");
    }
  }
  std::vector<XML_Node*> ccc;
  xctml->getChildren("simulation", ccc);
  int sz = ccc.size(); 
  if (num < 0 || num >= sz) {
    throw CanteraError("getSimulNum", "out of bounds");
  }
  return ccc[num];
}
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
  printf("  xmlSolnDiff [-h] [-a atol] [-r rtol] File1.xml File2.xml\n");
  printf("\t\n");
  printf("\tCompares two solution files in Cantera Solution XML format.\n");
  printf("\tThe comparison is done using a weighted norm basis.\n");
  printf("\t\n");
  printf("\tThe two files should be basically equal. However, File1.xml is\n");
  printf("\ttaken as the reference file that has precedence, when there is\n");
  printf("\tsomething to be decided upon.\n");
  printf("\t\n");
  printf("\t Arguments:\n");
  printf("\t  -h      = Usage info\n");
  printf("\t  -a atol = Set absolute tolerance parameter - default = 1.0E-9\n");
  printf("\t  -r rtol = Set relative tolerance parameter - default = 1.0E-3\n");
  printf("\t\n");
  printf("\t Shell Return Values:\n");
  printf("\t   1 = Comparison was successful\n");
  printf("\t   0 = One or more nodal values failed the comparison\n");
  printf("\t  -1 = Apples to oranges, the files can not even be compared against\n");
  printf("\t       one another.\n");
  printf("\t\n");
}
//====================================================================================================================
int main(int argc, char *argv[])
{
  int opt_let;
  char  *fileName1=NULL, *fileName2=NULL;
  XML_Node *fp1 = 0;
  XML_Node *fp2 = 0;
  int    testPassed = 1;
  double atol_arg = 0.0, rtol_arg = 0.0;
    int id = 0;
  int id2 = 0;
  char *ggg = 0;
  char *rrr = 0;
  int iTotal = 0;
  /*
   * Interpret command line arguments
   */
  /* Loop over each command line option */
  while((opt_let = getopt(argc, argv, "ha:r:")) != EOF) {

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
  printf("xmlSolnDiff: XML Soln File comparison utility program\n");
  printf("             Version $Revision: 5 $\n");
  printf("             Harry K. Moffat Div. 1516 Sandia National Labs\n");
  printf("           \n");
  printf("             First  Cantera XML Solution File = %s\n", fileName1);
  printf("             Second Cantera XML Solution File = %s\n", fileName2); 
  printf("\n");
  printf("             Absolute tol = %g\n", gatol);
  printf("             Relative tol = %g\n", grtol);
  printf("----------------------------------------------------------------------\n");
  printf("\n");  
  
  /*
   *  Open up the two xml Files #1 and #2
   */

  if (!(fp1 = readXML(fileName1 ))) {
    fprintf(stderr,"Error opening up file1, %s\n", fileName1);
    exit(-1);
  }
  if (!(fp2 = readXML(fileName2))) {
    fprintf(stderr, "Error opening up file2, %s\n", fileName2);
    exit(-1);
  }
  int numSim1 = countSimulations(fp1);
  int numSim2 = countSimulations(fp2);
  if (numSim1 <= 0) {
    printf("Number of Simulation Entries in file 1 is zero\n");
    testPassed = -1;
    exit(-1);
  }
  if (numSim1 != numSim2) {
    printf("Number of Simulation Entries differ: %d %d\n", numSim1, numSim2);
    testPassed = -1;
    exit(-1);
  }

  if (testPassed == -1) {
      iTotal = 1;
  }
  return iTotal;

}
//======================================================================================================================
