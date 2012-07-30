#include "KernelFunctions.h"

namespace PLMD {

Kernel* KernelRegister::create( const std::string& type, const KernelOptions& ko, const bool& needderivs ){
  Kernel* kernel;
  if( type=="uniform" ){
      kernel=dynamic_cast<Kernel*>( new UniformKernel(ko) ); 
  } else {
      return NULL;
  }
  if( needderivs && !kernel->hasderivatives ) return NULL;
  return kernel;
}


KernelOptions::KernelOptions( const std::vector<double>& at, const std::vector<double>& sig, const double& w, const bool& norm ):
pos(at), 
width(sig), 
height(w),
normalize(norm)
{
}

Kernel::Kernel( const KernelOptions& ko ):
center(ko.pos),
width(ko.width),
height(ko.height)
{
  if( center.size()==width.size() ) diagonal=true;
  else diagonal=false;
}

std::vector<unsigned> Kernel::getSupport( const std::vector<double>& dx ){
  plumed_assert( ndim()==dx.size() );
  std::vector<unsigned> support( dx.size() );
  for(unsigned i=0;i<dx.size();++i) support[i]=static_cast<unsigned>(ceil( getCutoff(width[i])/dx[i] ));
  return support;
}

double Kernel::getDeterminant() const {
  double vol;
  if(diagonal){
     vol=1; for(unsigned i=0;i<width.size();++i) vol*=width[i];
  } else {

  }
  return vol;
}

double Kernel::evaluate( const std::vector<bool>& pbc, const std::vector<double>& range, const std::vector<double>& pos ){
  plumed_assert( pbc.size()==ndim() && range.size()==ndim() && pos.size()==ndim() );

  double r2=0;
  if(diagonal){ 
     double s;
     for(unsigned i=0;i<ndim();++i){
         if( pbc[i] ){
             s=( pos[i] - center[i] )/range[i]; s=Tools::pbc(s); 
             s=( s*range[i] ) / width[i];
             r2+=s*s;
         } else {
             s=(pos[i] - center[i] )/ width[i];
             r2+=s*s;
         }
     }
  } else {

  }
  return getValue( sqrt(r2) );
}

UniformKernel::UniformKernel( const KernelOptions& ko ):
Kernel(ko)
{
  hasderivatives=false;
  if(ko.normalize){
    double vol;
    if( ko.pos.size()%2==1 ){
        double dfact=1;
        for(unsigned i=1;i<ko.pos.size();i+=2) dfact*=static_cast<double>(i);
        vol=( pow( pi, (ko.pos.size()-1)/2 ) ) * ( pow( 2., (ko.pos.size()+1)/2 ) ) / dfact;
    } else {
        double fact=1.;
        for(unsigned i=1;i<ko.pos.size()/2;++i) fact*=static_cast<double>(i);
        vol=pow( pi,ko.pos.size()/2 ) / fact;
    }
    vol*=getDeterminant();
    setHeight( ko.height/vol );
  }
}

}
