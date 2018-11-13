#include <stdio.h>


#include "cantera/base/logger.h"
#include "cantera/thermo.h"
#include "cantera/thermo/MineralEQ3.h"


#include <fstream>
#include <cstring>

using namespace Zuzax;
using namespace std;

void printUsage() {
    cout << "usage: findRefState [-h]  thermoFile.xml  " <<  endl;
    cout << "      " << endl;
}

int main(int argc, char** argv)
{
    int retn = 0;
    string commandFile = "xxx.xml";
    try {
        if (argc > 1) {
            std::string tok;
            for (int j = 1; j < argc; j++) {
                tok = string(argv[j]);
                if (tok[0] == '-') {
                    int nopt = static_cast<int>(tok.size());
                    for (int n = 1; n < nopt; n++) {
                        if (tok[n] == 'h') {
                            printUsage();
                            exit(1);
                        } else if (tok[n] == 'd') {
                            if (j < (argc - 1)) {
                                string tokla = string(argv[j+1]);
                                if (strlen(tokla.c_str()) > 0) {
                                    n = nopt - 1;
                                    j += 1;


                                }
                            }
                        } else {
                            printUsage();
                            exit(1);
                        }
                    }
                } else if (commandFile == "xxx.xml") {
                    commandFile = tok;
                } else {
                    printUsage();
                    exit(1);
                }
            }
        }

        string iFile = "O_Element.xml";
        iFile = commandFile;

        Zuzax::thermo_t_double* tp_ptr = Zuzax::newPhase(iFile, "");

        tp_ptr->setState_TP(298.15, Zuzax::OneBar);

        size_t nsp = tp_ptr->nSpecies();

        std::vector<double> S0_298(nsp, 0.0);
        std::vector<double> H298(nsp, 0.0);

        tp_ptr->getEntropy_R(S0_298.data());
        tp_ptr->getEnthalpy_RT(H298.data());

        for (double& val : S0_298) {
            val *= GasConstant / 1.0E3;
        }
        for (double& val : H298) {
            val *= GasConstant * 298.15 / 1.0E6;
        }

        printf("                 Species         S0(298)                          H298 \n");
        printf("                                (Joules/gmol/K)                   (kJ/gmol) \n");
        for (size_t k = 0; k < nsp; k++) {
            std::string ss = tp_ptr->speciesName(k);
            printf("  %24s  % 24.15E  % 24.15E \n", ss.c_str(), S0_298[k],
                   H298[k]);
        }
         printf("    % 24.15E   \n", S0_298[0]/2.0);

        return retn;

    } catch (Zuzax::ZuzaxError) {

        showErrors();
        appdelete();
        return -1;
    }
}

