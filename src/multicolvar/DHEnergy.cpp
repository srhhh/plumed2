/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "MultiColvar.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <iostream>
#include <string>

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC COLVAR DHENERGY2
/*
Calculate Debye-Huckel interaction energy among GROUPA and GROUPB. The two groups should
be disjointed.

\par Examples
\verbatim
# this is printing the electrostatic interaction between two groups of atoms
dh: DHEN GROUPA=1-10 GROUPB=11-20 EPS=80.0 I=0.1 TEMP=300.0
PRINT ARG=dh
\endverbatim
(see also \ref PRINT)

*/
//+ENDPLUMEDOC

class DHEnergy : public MultiColvar {
private:
  double k; // Inverse Debye screening length
  double constant;
  double epsilon;

public:
  static void registerKeywords( Keywords& keys );
  DHEnergy(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos );
/// Returns the number of coordinates of the field
  unsigned getNumberOfFieldDerivatives();
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(DHEnergy,"DHENERGY2")

void DHEnergy::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords(keys);
  ActionWithVessel::autoParallelize( keys );
  keys.use("GROUPA"); keys.use("GROUPB");
  keys.add("compulsory","I","1.0","Ionic strength (M)");
  keys.add("compulsory","TEMP","300.0","Simulation temperature (K)");
  keys.add("compulsory","EPSILON","80.0","Dielectric constant of solvent");
}

DHEnergy::DHEnergy(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao),
k(0.0),
constant(0.0)
{
  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // The components should just be summed
  addVessel("SUM","");  
  // And setup the ActionWithVessel
  requestDistribution();
  double I,T;
  parse("I",I);
  parse("TEMP",T);
  parse("EPSILON",epsilon);
  checkRead();

  plumed_assert(!plumed.getAtoms().usingNaturalUnits());
  atoms.getUnits().getLength();
  constant=138.935458111/atoms.getUnits().getEnergy();
  k=sqrt(I/(epsilon*T))*502.903741125;

  log<<"  with solvent dielectric constant "<<epsilon<<"\n";
  log<<"  at temperature "<<T<<" K\n";
  log<<"  at ionic strength "<<I<< "M\n";
  log<<"  Bibliography "<<plumed.cite("Trang, Carloni, Varani and Bussi, submitted (2013)")<<"\n";
}

unsigned DHEnergy::getNumberOfFieldDerivatives(){
  plumed_massert(0,"should not be in this routine");
  return 0;
}

double DHEnergy::compute( const unsigned& j, const std::vector<Vector>& pos ){
   Vector distance; distance=getSeparation( pos[0], pos[1] );

   if( getAbsoluteIndex(0)==getAbsoluteIndex(1) ) return 0.0;

   double value=distance.modulo();
   double invdistance=1.0/value;
   double tmp=exp(-k*value)*invdistance*constant*getCharge(0)*getCharge(1)/epsilon;
   double dtmp=-(k+invdistance)*tmp;

   // And finish the calculation
   addAtomsDerivatives( 0,-invdistance*dtmp*distance );
   addAtomsDerivatives( 1, invdistance*dtmp*distance );
   addBoxDerivatives( -invdistance*dtmp*Tensor(distance,distance) );
   return tmp;
}

}
}
