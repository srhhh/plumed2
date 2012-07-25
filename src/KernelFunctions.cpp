#include "KernelFunctions.h"

namespace PLMD {

std::vector<unsigned> Kernel::getSupport( const std::vector<double>& dx ){
  plumed_assert( ndim()==dx.size() );
  std::vector<unsigned> support( dx.size() );
  for(unsigned i=0;i<dx.size();++i) support[i]=static_cast<unsigned>(ceil( getCutoff(width[i])/dx[i] ));
  return support;
}

double Kernel::evaluate( const std::vector<bool>& pbc, const std::vector<double>& range, const std::vector<double>& pos ){
  plumed_assert( pbc.size()==ndim() && range.size()==ndim() && pos.size()==ndim() );
  std::vector<double> dx( pos.size() ); double s;
  for(unsigned i=0;i<ndim();++i){
      if( pbc[i] ){
          s=( pos[i] - center[i] )/range[i];
          s=Tools::pbc(s); dx[i]=( s*range[i] ) / width[i];
      } else {
          dx[i]=(pos[i] - center[i] )/ width[i];
      }
  }
  return getValue( dx );
}

UniformKernel::UniformKernel( const std::vector<double>& at, const std::vector<double>& sig, const double& w ):
Kernel( at,sig,w )
{
  double vol=1; for(unsigned i=0;i<sig.size();++i) vol*=0.5*sig[i];
  setHeight( w/vol );
}

}
