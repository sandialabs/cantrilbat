#include "exp_BoundaryCondition.h"
#include "exp_DakotaInterface.h"

using namespace Cantera;


int main() {
  
  try {
    BCsteptable bc( "sampleBC.xml" );
  }
  catch (CanteraError) {
    showErrors(std::cout);
  }
  
}
