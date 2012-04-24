#include "MultiColvar.h"
#include "ActionRegister.h"
#include "SwitchingFunction.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC MCOLVAR CONTACTMAP
/**

*/
//+ENDPLUMEDOC

class ContactMap : public MultiColvar {
private:
  double rcut;
  std::vector<SwitchingFunction> sfs;
public:
  static void registerKeywords( Keywords& keys );
  ContactMap(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
/// Returns the number of coordinates of the field
  unsigned getNumberOfFieldDerivatives();
  bool isPeriodic(const unsigned nn){ return false; }
/// Setup the field
  void derivedFieldSetup( const double sigma ){}
/// Calculate the contribution to the field at the point thisp
  void calculateFieldContribution( const unsigned& j, const std::vector<double>& thisp, Value* tmpvalue, Value& tmpstress, std::vector<Value>& tmpder );
};

PLUMED_REGISTER_ACTION(ContactMap,"CONTACTMAP")

void ContactMap::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithDistribution::autoParallelize( keys );
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("numbered","SWITCH","The switching functions to use for each of the contacts in your map.  You can either specify a global switching function using SWITCH or one switching function for each contact"); keys.reset_style("SWITCH","compulsory"); 
  keys.use("ATOMS");
}

ContactMap::ContactMap(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao),
rcut(-1)
{
  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // Read the cutoff for the neighbour list
  if( isTimeForNeighborListUpdate() ){
      parse("NL_CUTOFF",rcut);
      if( rcut>0 ) log.printf("  ignoring distances greater than %lf in neighbor list\n",rcut);
  }
  // Functon is not periodic
  setNotPeriodic();

  // Read in switching functions
  std::string sw, errors; parse("SWITCH",sw);
  for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i) sfs.push_back( SwitchingFunction() );
  if(sw.length()>0){
     for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){
        sfs[i].set(sw,errors);
        if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
     }
  } else {
     for(unsigned i=0;i<getNumberOfFunctionsInDistribution();++i){
         std::string num, sw1; Tools::convert(i+1, num);
         if( !parseNumbered( "SWITCH", i+1, sw1 ) ) error("missing SWITCH" + num + " keyword");
         sfs[i].set(sw1,errors);
         if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
     }
  }

  // And check everything has been read in correctly
  checkRead();
  // Create the field
  addField("", new Field( "cvlist",1) );
  // And setup the ActionWithDistribution
  requestDistribution();
}

unsigned ContactMap::getNumberOfFieldDerivatives(){
  return getNumberOfFunctionsInDistribution();
}

double ContactMap::compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
   Vector distance;
   distance=getSeparation( pos[0], pos[1] );
   double value, dfunc;
   value=sfs[j].calculate( distance.modulo() , dfunc);

   // Check at neighbor list update time whether this distance is big
   if( isTimeForNeighborListUpdate() && rcut>0 ){
       if( value>rcut ){ stopCalculatingThisCV(); return 0.0; }
   }

   // And finish the calculation
   deriv[0]=-dfunc*distance;
   deriv[1]=dfunc*distance;
   virial=-dfunc*Tensor(distance,distance);
   return value;
}

void ContactMap::calculateFieldContribution( const unsigned& j, const std::vector<double>& thisp, Value* tmpvalue, Value& tmpstress, std::vector<Value>& tmpder ){
   if( static_cast<double>(j)==thisp[0] ) tmpstress.set( tmpvalue->get() );
   else tmpstress.set(0);
//   printf("ARSE HELLO %d %f %f %f \n",j,thisp[0],tmpvalue->get(), tmpstress.get() );
}

}
