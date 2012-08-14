#include "KernelFunctions.h"

namespace PLMD {

Kernel::Kernel( const std::vector<double>& at, const std::vector<double>& sig, const std::string& type, const double& w, const bool& norm ):
center(at),
width(sig)
{
  unsigned ncv=center.size();
  if( width.size()==ncv ) diagonal=true;
  else if( width.size()==(ncv*(ncv+1))/2 ) diagonal=false;
  else plumed_massert(0,"specified sigma is neither diagonal or full covariance matrix");

  // Setup the non linear function
  std::vector<std::string> params;
  nlfunc.set(type, params );

  if( norm ){
    double det;
    if(diagonal){
       det=1; for(unsigned i=0;i<width.size();++i) det*=width[i];
    } else {
       unsigned ncv=ndim(); 
       Matrix<double> mymatrix( getMatrix() ), myinv( ncv, ncv );
       Invert(mymatrix,myinv); double logd;
       logdet( myinv, logd );
       det=exp(logd);
    }
    height=w / nlfunc.getVolume(ncv, det );
  } else {
    height=w;
  }
}

Kernel::Kernel( const std::vector<std::string>& cv_names, PlumedIFile& ifile ){
  center.resize( cv_names.size() );
  // Read the position of the center
  for(unsigned i=0;i<cv_names.size();++i) ifile.scanField( cv_names[i], center[i] );
  // Read the covariance
  std::string sss; ifile.scanField("multivariate",sss);
  if(sss=="true") diagonal=false;
  else if(sss=="false") diagonal=true;
  else plumed_merror("cannot parse multivariate = " +sss );
  if( diagonal || cv_names.size()==1 ){
      width.resize( cv_names.size() );
      for(unsigned i=0;i<cv_names.size();++i) ifile.scanField("sigma_"+cv_names[i],width[i]);
  } else {
      unsigned ncv=cv_names.size();
      width.resize( (ncv*(ncv+1))/2 );
      Matrix<double> upper(ncv,ncv), lower(ncv,ncv);
      for (unsigned i=0;i<ncv;i++){
          for (unsigned j=0;j<ncv-i;j++){
              ifile.scanField("sigma_"+cv_names[j+i]+"_"+cv_names[j],lower(j+i,j));
              upper(j,j+i)=lower(j+i,j);
          }
      }
      Matrix<double> mymult(ncv,ncv), invmatrix(ncv,ncv);
      mult(lower,upper,mymult);
      // now invert and get the sigmas
      Invert(mymult,invmatrix);
      unsigned k=0;
      for (unsigned i=0;i<ncv;i++){
           for (unsigned j=i;j<ncv;j++){
                width[k]=invmatrix(i,j);
                k++;
           }
      }
  }
  // Read the height
  ifile.scanField("height",height);
  // Read the kernel type
  std::string ktype; ifile.scanField( "kerneltype", ktype );
  nlfunc.set( ktype, ifile );
}

std::vector<unsigned> Kernel::getSupport( const std::vector<double>& dx ) const {
  plumed_assert( ndim()==dx.size() );
  std::vector<unsigned> support( dx.size() );
  if(diagonal){
     for(unsigned i=0;i<dx.size();++i) support[i]=static_cast<unsigned>(ceil( nlfunc.getCutoff(width[i])/dx[i] ));
  } else {
     unsigned ncv=ndim(); 
     Matrix<double> mymatrix( getMatrix() ), myinv( ncv,ncv );
     Invert(mymatrix,myinv);
     Matrix<double> myautovec(ncv,ncv); std::vector<double> myautoval(ncv);  
     diagMat(myinv,myautoval,myautovec);
     for(unsigned i=0;i<dx.size();++i){
         double extent=fabs(sqrt(myautoval[0])*myautovec(i,0)); 
         support[i]=static_cast<unsigned>(ceil( nlfunc.getCutoff( extent )/dx[i] ));
     }
  }
  return support;
}

double Kernel::evaluate( const std::vector<Value*>& pos, std::vector<double>& derivatives, bool usederiv ) const {
  plumed_assert( pos.size()==ndim() && derivatives.size()==ndim() );
  if( usederiv ) plumed_assert( nlfunc.hasDerivatives() ); 

  double r2=0;
  if(diagonal){ 
     for(unsigned i=0;i<ndim();++i){
         derivatives[i]=pos[i]->difference( center[i] ) / width[i]; 
         r2+=derivatives[i]*derivatives[i];
     }
  } else {
     Matrix<double> mymatrix( getMatrix() ); 
     for(unsigned i=0;i<mymatrix.nrows();++i){
        double dp_i, dp_j; derivatives[i]=0;
        dp_i=pos[i]->difference( center[i] ); 
        for(unsigned j=0;j<mymatrix.ncols();++j){
          if(i==j) dp_j=dp_i;
          else dp_j=pos[j]->difference( center[j] );

          derivatives[i]+=mymatrix(i,j)*dp_j;
          r2+=dp_i*dp_j*mymatrix(i,j);
        }
     }
  }
  double kderiv, kval;
  
  if( !nlfunc.inputXSquared() ){
     double r=sqrt(r2);
     if(!usederiv){
         kval=height*nlfunc.calculate(r); 
     } else {
         kval=height*nlfunc.calculate( r, kderiv ); kderiv*=height / r;
     }
  } else {
     if(!usederiv){
         kval=height*nlfunc.calculate(r2);
     } else {
         kval=height*nlfunc.calculate( r2, kderiv ); kderiv*=height;
     }
  }
  if(usederiv){
     for(unsigned i=0;i<ndim();++i) derivatives[i]*=kderiv;
  }
  return kval;
}

void Kernel::print( const std::vector<std::string>& cv_names, PlumedOFile& ofile ) const {
  plumed_assert( cv_names.size()==ndim() );
  for(unsigned i=0;i<ndim();++i){ ofile.printField( cv_names[i],center[i] ); }   //fprintf(ofile, "%14.9f   ", center[i]);
  if(ndim()==1){
      ofile.printField( "multivariate","false" );
      ofile.printField( "sigma_" + cv_names[0], width[0]);
  } else if(diagonal){
      ofile.printField( "multivariate","false" );
      for(unsigned i=0;i<ndim();++i) ofile.printField( "sigma_" + cv_names[i], width[i] );     // { fprintf(ofile, "%14.9f   ", width[i]); }
  } else{
      ofile.printField( "multivariate","true");
      Matrix<double> mymatrix( getMatrix() );
      // invert the matrix
      Matrix<double> invmatrix( ndim(),ndim() );
      Invert(mymatrix,invmatrix);
      // enforce symmetry
      for(unsigned i=0;i<ndim();i++){
          for(unsigned j=i;j<ndim();j++){
              invmatrix(i,j)=invmatrix(j,i);
          }
      }
      Matrix<double> lower( ndim() ,ndim() );
      cholesky(invmatrix,lower); // now this , in band form , is similar to the sigmas
      // loop in band form 
      unsigned k=0;
      for(unsigned i=0;i<ndim();i++){
          for(unsigned j=0;j<ndim()-i;j++){
              ofile.printField("sigma_" + cv_names[j+i] + "_" + cv_names[j], lower(j+i,j));
              k++;
          }
      }
  }
  ofile.printField("height",height);
  nlfunc.printParameters(ofile);
}

}
