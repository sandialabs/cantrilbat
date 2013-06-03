#include "exp_BoundaryCondition.h"
#include "exp_DakotaInterface.h"

using namespace Cantera;


int main() {

  DakotaInterface di;
  std::cout << "has " << di.hasParams() <<  " parameters" << std::endl;
  if ( di.hasParam( "current" ) )
    std::cout << "Found current "<< std::endl << std::endl;

  di.writeParameters();

}
