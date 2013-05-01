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
#include "SingleDomainRMSD.h"
#include "DRMSD.h"

namespace PLMD {

SingleDomainRMSD::SingleDomainRMSD( const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration(ro),
ReferenceAtoms(ro)
{
}

void SingleDomainRMSD::setReference( const std::vector<Vector>& conf, const std::vector<double> align_in, const std::vector<double> displace_in ){
  reference_atoms.resize( conf.size() ); align.resize( conf.size() ); 
  displace.resize( conf.size() ); der_index.resize( conf.size() );
  for(unsigned i=0;i<conf.size();++i){
     reference_atoms[i]=conf[i]; align[i]=align_in[i]; 
     displace[i]=displace_in[i]; der_index[i]=i;
  }
  setNumberOfAtoms( conf.size() ); setNumberOfArguments( 0 );
}

double SingleDomainRMSD::calculate( const std::vector<Vector>& pos, const Pbc& pbc,  const bool& squared ){
  clearDerivatives();
  return calc( pos, pbc, squared );
}

double SingleDomainRMSD::calc( const std::vector<Vector>& pos, const Pbc& pbc, const std::vector<Value*>& vals, const std::vector<double>& arg, const bool& squared ){
  plumed_dbg_assert( vals.size()==0 && pos.size()==getNumberOfAtoms() && arg.size()==0 );
  return calc( pos, pbc, squared );
}

}
