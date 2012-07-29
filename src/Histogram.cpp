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
#include "Analysis.h"
#include "ActionRegister.h"
#include "Grid.h"
#include "KernelFunctions.h"

namespace PLMD{

//+PLUMEDOC ANALYSIS HISTOGRAM
/* 
Calculate a histogram as a function of a few CVs

\par Examples

*/
//+ENDPLUMEDOC

class CVHistogram : public Analysis {
private:
  Grid* gg;
  FILE* gridfile;
  double normalization;
  std::vector<double> point, bw;
  std::string kerneltype; 
  unsigned wgridstride;
public:
  static void registerKeywords( Keywords& keys );
  CVHistogram(const ActionOptions&ao);
  ~CVHistogram();
  void performAnalysis();
};

PLUMED_REGISTER_ACTION(CVHistogram,"HISTOGRAM")

void CVHistogram::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("compulsory","GRID_BIN","the number of bins for the grid");
  keys.add("compulsory","KERNEL","the kernel function you are using");
  keys.add("compulsory","BANDWIDTH","the bandwdith for kernel density estimation");
  keys.add("compulsory","GRID_WSTRIDE","write the histogram to a file every N steps");
  keys.add("compulsory","GRID_WFILE","histogram","the file on which to write the grid");
}

CVHistogram::CVHistogram(const ActionOptions&ao):
PLUMED_ANALYSIS_INIT(ao),
normalization(0),
point(getNumberOfArguments()),
bw(getNumberOfArguments())
{
  // Read stuff for Grid
  std::vector<double> gmin(getNumberOfArguments());
  parseVector("GRID_MIN",gmin);
  std::vector<double> gmax(getNumberOfArguments());
  parseVector("GRID_MAX",gmax);
  std::vector<unsigned> gbin(getNumberOfArguments());
  parseVector("GRID_BIN",gbin);
  parse("GRID_WSTRIDE",wgridstride);
  if( wgridstride%getStride()!=0 ) error("frequency of grid writing must be a multiple of the frequency for data collection");
  std::string gridfname; parse("GRID_WFILE",gridfname); 

  // Read stuff for window functions
  parseVector("BANDWIDTH",bw);
  // Read the type of kernel we are using
  parse("KERNEL",kerneltype);

  std::vector<double> point, bw; 
  Kernel* kernel=KernelRegister::create( kerneltype, KernelOptions(point, bw, 1.0, true), false );
  if(!kernel) error("not a valid kernel function " + kerneltype );
  delete kernel; 
  checkRead();

  log.printf("  Grid min");
  for(unsigned i=0;i<gmin.size();++i) log.printf(" %f",gmin[i]);
  log.printf("\n");
  log.printf("  Grid max");
  for(unsigned i=0;i<gmax.size();++i) log.printf(" %f",gmax[i]);
  log.printf("\n");
  log.printf("  Grid bin");
  for(unsigned i=0;i<gbin.size();++i) log.printf(" %d",gbin[i]);
  log.printf("\n");

  // Get pbc stuff for grid
  std::vector<bool> pbc;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    pbc.push_back(getPntrToArgument(i)->isPeriodic());
    if(pbc[i]){
     double dmin,dmax;
     getPntrToArgument(i)->getDomain(dmin,dmax);
     gmin[i]=dmin; gmax[i]=dmax;
    }
  }

  // Create the grid
  gg=new Grid(gmin,gmax,gbin,pbc,false,false);
  // And a file to write it on
  gridfile=fopen(gridfname.c_str(),"w");
}

CVHistogram::~CVHistogram(){
  delete gg;
  if(gridfile) fclose(gridfile); 
}

void CVHistogram::performAnalysis(){
  double weight; getDataPoint( 0, point, weight );
  normalization+=weight;
  Kernel* kernel=KernelRegister::create( kerneltype, KernelOptions(point, bw, weight, true), false );
  gg->addKernel( kernel );
  delete kernel;

//  if( kernel==uniform ){
//      UniformKernel* unif=new UniformKernel( KernelOptions(point, bw, weight,true) );
//      Kernel* kern=dynamic_cast<Kernel*>( unif );
//      gg->addKernel( kern );
//  } else if( kernel==triangular ){
////      TriangularKernel triangle( point, bw, weight );
////      gg->add( triangle );
//  }

// dump grid on file
  if(wgridstride>0 && getStep()%wgridstride==0){
     gg->writeToFile(gridfile,normalization);
  }  
}

}
