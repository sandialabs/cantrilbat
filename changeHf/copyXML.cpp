/**
 *  @file copyXML.cpp
 *
 */

//  Example 
//
//  This does a complete read and then copy of an
//  XML file.


#include <cstdio>
#include <fstream>

using namespace std;
using namespace Cantera;

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
static void printUsage()
{
    printf("copyXML infile.xml outfile.xml -h\n");
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

int main(int argc, char** argv) {
    string infile = "";
    string outfile = "";
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
	} else if (infile == "") {
	  infile = tok;
	} else if (outfile == "") {
	  outfile = tok;
	}
	else {
	  printf("copyXML ERROR: too many tokens on input line\n");
	  printUsage();
	  exit(1);
	}
      }
    }

    if (infile.size() <= 0) {
      printf("copyXML ERROR: must specify input XML file name\n");
      printUsage();
      exit(-1);
    }
    if (outfile.size() <= 0) {
      printf("copyXML ERROR: must specify output XML file name\n");
      printUsage();
      exit(-1);
    }
 
    XML_Node *xc = new XML_Node();
    const char *fn = infile.c_str();
    string path = findInputFile(fn);
    ifstream fin(path.c_str());
    if (!fin) {
      throw ZuzaxError("copyXML","could not open "
			 +path+" for reading.");
    }
    /*
     * Make a complete copy of the xml file
     */
    xc->build(fin);
    fin.close();
    XML_Node *xd = new XML_Node();
    xc->copy(xd);
 
    ofstream tout;
    tout.open(outfile.c_str());
    xc->write(tout);
    tout.close();

    delete xc;
    delete xd;
    appdelete();

    return 0;
}
/***********************************************************/
