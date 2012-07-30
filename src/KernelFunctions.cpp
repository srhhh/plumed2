#include "KernelFunctions.h"

namespace PLMD {

Kernel* KernelRegister::create( const std::string& type, const KernelOptions& ko, const bool& needderivs ){
  Kernel* kernel;
  if( type=="uniform" ){
      kernel=dynamic_cast<Kernel*>( new UniformKernel(ko) ); 
  } else if( type=="gaussian"){
      kernel=dynamic_cast<Kernel*>( new GaussianKernel(ko) );
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
height(ko.height),
is_function_of_r2(false)
{
  unsigned ncv=center.size();
  if( width.size()==ncv ) diagonal=true;
  else if( width.size()==(ncv*(ncv+1))/2 ) diagonal=false;
  else plumed_massert(0,"specified sigma is neither diagonal or full covariance matrix");
}

std::vector<unsigned> Kernel::getSupport( const std::vector<double>& dx ){
  plumed_assert( ndim()==dx.size() );
  std::vector<unsigned> support( dx.size() );
  if(diagonal){
     for(unsigned i=0;i<dx.size();++i) support[i]=static_cast<unsigned>(ceil( getCutoff(width[i])/dx[i] ));
  } else {
     unsigned ncv=ndim(); 
     Matrix<double> mymatrix( getMatrix() ), myinv( ncv,ncv );
     Invert(mymatrix,myinv);
     Matrix<double> myautovec(ncv,ncv); std::vector<double> myautoval(ncv);  
     diagMat(myinv,myautoval,myautovec);
     for(unsigned i=0;i<dx.size();++i){
         double extent=fabs(sqrt(myautoval[0])*myautovec(i,0)); 
         support[i]=static_cast<unsigned>(ceil( getCutoff( extent )/dx[i] ));
     }
  }
  return support;
}

double Kernel::getDeterminant() const {
  double vol;
  if(diagonal){
     vol=1; for(unsigned i=0;i<width.size();++i) vol*=width[i];
  } else {
     unsigned ncv=ndim(); 
     Matrix<double> mymatrix( getMatrix() ), myinv( ncv, ncv );
     Invert(mymatrix,myinv); double logd;
     logdet( myinv, logd );
     vol=exp(logd);
  }
  return vol;
}

double Kernel::evaluate( const std::vector<Value>& pos, std::vector<double>& derivatives, bool usederiv ){
  plumed_assert( pos.size()==ndim() && derivatives.size()==ndim() );
  if( usederiv ) plumed_assert( hasderivatives ); 

  double r2=0;
  if(diagonal){ 
     for(unsigned i=0;i<ndim();++i){
         derivatives[i]=pos[i].difference( center[i] ) / width[i]; 
         r2+=derivatives[i]*derivatives[i];
     }
  } else {
     Matrix<double> mymatrix( getMatrix() ); 
     for(unsigned i=0;i<mymatrix.nrows();++i){
        double dp_i, dp_j; derivatives[i]=0;
        dp_i=pos[i].difference( center[i] ); 
        for(unsigned j=0;j<mymatrix.ncols();++j){
          if(i==j) dp_j=dp_i;
          else dp_j=pos[j].difference( center[j] );

          derivatives[i]+=mymatrix(i,j)*dp_j;
          r2+=dp_i*dp_j*mymatrix(i,j);
        }
     }
  }
  double kderiv, kval;
  
  if( !is_function_of_r2 ){
     double r=sqrt(r2);
     kval=getValue( r, kderiv ); kderiv/=r;
  } else {
     kval=getValue( r2, kderiv );
  }
  for(unsigned i=0;i<ndim();++i) derivatives[i]*=kderiv;
  return kval;
}

UniformKernel::UniformKernel( const KernelOptions& ko ):
Kernel(ko)
{
  hasderivatives=false; is_function_of_r2=false;
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

GaussianKernel::GaussianKernel( const KernelOptions& ko ):
Kernel(ko),
DP2CUTOFF(6.25)
{
  hasderivatives=true; is_function_of_r2=true;
  if(ko.normalize){
     double vol;
     vol=( pow( 2*pi, 0.5*ko.pos.size() ) * pow( getDeterminant(), 0.5 ) );
     setHeight( ko.height/vol);
  }
}

}
