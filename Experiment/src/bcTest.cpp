#include "exp_BoundaryCondition.h"
#include "exp_DakotaInterface.h"

using namespace Zuzax;


int main() {
  
  try {
    BCsteptable bc( "sampleBC.xml" );
  }
  catch (ZuzaxError) {
    showErrors(std::cout);
  }
} 
