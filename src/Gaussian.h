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

#include <vector>
using namespace std;

/// a simple structure to conveniently store Gaussian from Metadynamics 
struct Gaussian {
 vector<double> center;
 vector<double> sigma;
 double height;
 bool   multivariate; // this is required to discriminate the one dimensional case 
 vector<double> invsigma;
 Gaussian(const vector<double> & center,const vector<double> & sigma,double height, bool multivariate ):
   center(center),sigma(sigma),height(height),multivariate(multivariate),invsigma(sigma){
     for(unsigned i=0;i<invsigma.size();++i)abs(invsigma[i])>1.e-20?invsigma[i]=1.0/invsigma[i]:0.; // to avoid troubles from zero element in flexible hills
   }
};
 
