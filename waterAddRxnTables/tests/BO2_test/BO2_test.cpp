/*
 *
 */

#include "zuzax/IdealGasMix.h"
#include "zuzax/equilibrium.h"
#include "zuzax/multiphase/MP_EquilStatic.h"

#include "zuzax/thermo.h"
#include "zuzax/equil/vcs_internal.h"
#include "zuzax/base/logger.h"
#include "zuzax/thermo/HMWSoln.h"
#include "zuzax/runners/pH_Speciation.h"

#include <cstring>

using namespace Zuzax;
using namespace std;

//=================================================================================================================================
void printUsage() {
    cout << "usage: nacl_dessicate [-h] [-help_cmdfile] [-d #] "
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << endl;
}
//=================================================================================================================================
void checkReturn(size_t ii, std::string ss = "" )
{
   if (ii == npos) {
        throw ZuzaxError("checkReturn()", "index not found" + ss);
   }
}
//=================================================================================================================================

//=================================================================================================================================
class my_pH_Speciation  : public pH_Speciation
{
public:

    //! Default constructor
    my_pH_Speciation() :
         pH_Speciation()
    {

        if (ip) {
            delete ip;
        }
        ip = new my_inputParams();
    }


    //! Default destructor
    virtual ~my_pH_Speciation()
    {
    }

    class my_inputParams : public inputParams
    {
    public:
        //! constructor
        /*!
         *  This is where we can override the base inputParams conditions
         */
        my_inputParams() :
          inputParams()
        {
            //  inputFile = "HMW.xml";
            //  PLinputFile {"---"};
            //! Type of initial brine composition
            /*!
             *    composition types
             *        0                   String version of Molarities of solutes
             *        1                   String version of Molalities of solutes
             *        2                   String version of Moles of all Species in the mixtures
             */ 
            brineCompositionType = 0;
            initialBrineComposition = "Na+:0.0101 Cl-:0.011 OH-:8.0E-05 H+:2.0E-09 UO2++:1.0E-4 U++++:1.0E-30 UO2CO3(aq):1.0E-10 UO2(CO3)2--:1.0E-10"
                                      " UO2(CO3)3-4:1.0E-10 HCO3-:7.0E-05 O2(aq):1.0E-10 H2(aq):1.0E-10";

            // stringPhaseMoles = "  " ;
            //spInc_Molality = "specialElement HCO3-";
            spInc_Molality = "all";

        }

        my_inputParams(const struct inputParams& rs) :
           inputParams()     
        {
            operator=(rs);
        }

        my_inputParams& operator=(const my_inputParams& rs)
        {
            if (this == &rs) {
                return *this;
            }
            inputParams::operator=(rs);
        
            return *this;
        }

        virtual ~my_inputParams()
        {
        }
    };


};


//=================================================================================================================================
int main(int argc, char **argv) {
  try {
    //int solver = 2;
    
    vcs_nonideal::vcs_timing_print_lvl = 0;
    Zuzax::s_VCS_Write_CSV_Report = 1;
    //int printLvl = 2;
    //int estimateEquil = 0;
  /*
   * Process the command line arguments
   */
  if (argc > 1) {
    string tok;
    for (int j = 1; j < argc; j++) {
      tok = string(argv[j]);
      if (tok[0] == '-') {
	int nopt = static_cast<int>(tok.size());
	for (int n = 1; n < nopt; n++) {
	  if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
	  } else if (tok[n] == 'h') {
	    printUsage();
	    exit(1);
	  } else if (tok[n] == 'd') {
	    //printLvl = 2;
	    int lvl = 2;
	    if (j < (argc - 1)) {
	      string tokla = string(argv[j+1]);
	      if (strlen(tokla.c_str()) > 0) {
		lvl = atoi(tokla.c_str());
		n = nopt - 1;
		j += 1;
		if (lvl >= 0 && lvl <= 1000) {
		  //if (lvl == 0) printLvl = 0;
		  //else          printLvl = lvl;
		}
	      }
	    }
	  } else {
	    printUsage();
	    exit(1);
	  }
	}
      } else {
	printUsage();
	exit(1);
      }
    }
  }

    my_pH_Speciation eRunner;


    //! Set up the reactor network
    eRunner.setupPhaseList();

    //! Initialize the network
    eRunner.initialize();


    eRunner.runCalculation();


    eRunner.writeResults();
    
  }
  catch (ZuzaxError& zee) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
