#include "NonLinearFunctions.h"

using namespace PLMD;

NonLinearFunction::NonLinearFunction():
setup(false)
{
}

std::string NonLinearFunction::printLatexFunctions( const std::string& x, const bool needderiv ){
  std::ostringstream ostr;
  if(!needderiv){
  
  }   
  ostr<<"\\f$s(r)=\\frac{ 1 - \\left( "<<x<<" \\right)^{n} }{ 1 - \\left( "<<x<<" \\right)^{m} } \\f$, ";
  ostr<<"\\f$s(r)=\\exp\\left(-"<<x<<"\\right)\\f$, \\f$s(r)=\\exp\\left(-"<<x<<"\\right)\\f$, ";
  ostr<<"or using \\f$s(r)=1-\\left|"<<x<<"\\right|\\f$.";
  return ostr.str();
}

std::string NonLinearFunction::printInput( const std::string& aparams, const bool needderiv ){
  std::ostringstream ostr;
  if(!needderiv){

  } else {
     ostr<<"The first of these options is specified using the syntax ";
  }
  ostr<<"{RATIONAL "<<aparams<<" NN=\\f$n\\f$ MM=\\f$m\\f$} and if ";
  ostr<<"the NN and MM keywords are missing they are assumed equal to 6 and 12 respectively. The second form is specified using ";
  ostr<<"{EXP "<<aparams<<"}. The third form is specified using {GAUSSIAN "<<aparams<<"}. The fourth form is specified using ";
  ostr<<"{TRIANGULAR "<<aparams<<"}.";
  return ostr.str();
}

void NonLinearFunction::set( const std::string& name, std::vector<std::string>& data ){
  setup=true;

  if(name=="RATIONAL" || name=="rational"){
      type=rational;
      nn=6; Tools::parse(data,"NN",nn);
      mm=12; Tools::parse(data,"MM",mm);
  } else if(name=="EXP" || name=="exp"){
      type=exponential;
  } else if(name=="GAUSSIAN" || name=="gaussian"){
      type=gaussian;
  } else if(name=="STEP" || name=="step"){
      type=step;
  } else if(name=="TRIANGULAR" || name=="triangular"){
      type=triangular;
  } else {
      plumed_merror( "Not a valid function type " + name );
  } 
}

void NonLinearFunction::set( const std::string& name, PlumedIFile& ifile ){
  setup=true;
  
  if(name=="RATIONAL" || name=="rational"){
      type=rational;
      ifile.scanField("NN",nn);  // nn=6; Tools::parse(data,"NN",nn);
      ifile.scanField("MM",mm);  // mm=12; Tools::parse(data,"MM",mm);
  } else if(name=="EXP" || name=="exp"){
      type=exponential;
  } else if(name=="GAUSSIAN" || name=="gaussian"){
      type=gaussian;
  } else if(name=="STEP" || name=="step"){
      type=step;
  } else if(name=="TRIANGULAR" || name=="triangular"){
      type=triangular;
  } else {
      plumed_merror( "Not a valid function type " + name );
  }
}

double NonLinearFunction::calculate( const double& x ) const {
  plumed_massert(setup,"non linear function has not been set"); 
  if(type==rational){
      if(x>(1.-100.0*epsilon) && x<(1+100.0*epsilon)){
         return nn/mm;
      }else{
         double rNdist=x;
         double rMdist=x;
    // this is a naive optimization
    // we probably have to implement some generic, fast pow(double,int)
         if(nn>2) for(int i=0;i<nn-2;i++) rNdist*=x;
         else rNdist = pow(x, nn-1);
         if(mm>2) for(int i=0;i<mm-2;i++) rMdist*=x;
         else rMdist = pow(x, mm-1);
         double num = 1.-rNdist*x;
         double iden = 1./(1.-rMdist*x);
         return num*iden;
      }
  } else if(type==exponential){
      return exp(-x);
  } else if(type==gaussian){
      return exp(-0.5*x);
  } else if(type==triangular){
      if( x<1.0 ){
         return  1. - fabs(x);
      }
      return 0;
  } else if(type==step){
      if(x<1.0) return 1;
      return 0;
  } else {
     plumed_merror("not a valid function type");
  }
  return 0;
}

double NonLinearFunction::calculate( const double& x, double& dx ) const {
  plumed_massert(setup,"non linear function has not been set");

  if(type==rational){
      if(x>(1.-100.0*epsilon) && x<(1+100.0*epsilon)){
         dx=0.5*nn*(nn-mm)/mm;
         return nn/mm;
      }else{
         double rNdist=x;
         double rMdist=x;
    // this is a naive optimization
    // we probably have to implement some generic, fast pow(double,int)
         if(nn>2) for(int i=0;i<nn-2;i++) rNdist*=x;
         else rNdist = pow(x, nn-1);
         if(mm>2) for(int i=0;i<mm-2;i++) rMdist*=x;
         else rMdist = pow(x, mm-1);
         double num = 1.-rNdist*x;
         double iden = 1./(1.-rMdist*x);
         double func = num*iden;
         dx = ((func*(iden*mm)*rMdist)-(nn*rNdist*iden));
         return func;
      }
  } else if(type==exponential){
      double ans=exp(-x);
      dx=-ans;
      return ans;
  } else if(type==gaussian){
      double ans=exp(-0.5*x);
      dx=-ans;        /// N.B. Here we have the derivative with respect to the square
      return ans;
  } else if(type==triangular){
      if( x<1.0 ){
         if(x==0) dx=0;
         else if(x>0) dx=-1; 
         else dx=1;
         return  1. - fabs(x); 
      }
      dx=0.0;
      return 0;
  } else {
      plumed_merror("derivatives unavailable for function type");
  } 
  dx=0;
  return 0;    
}

std::string NonLinearFunction::getName() const {
  plumed_massert(setup,"non linear function has not been set");  

  if(type==rational) return "rational";
  else if(type==exponential) return "exponential";
  else if(type==gaussian) return "gaussian";
  else if(type==triangular) return "triangular";
  else if(type==step) return "step";
  return "";
}

std::string NonLinearFunction::writeParameters() const {
  plumed_massert(setup,"non linear function has not been set");

  if(type==rational){
     std::ostringstream ostr; ostr<<"nn="<<nn<<" mm="<<mm;
     return ostr.str(); 
  }
  return "";
}

void NonLinearFunction::printParameters( PlumedOFile& ofile ) const {
  plumed_massert(setup,"non linear function has not been set");

  if(type==rational){
     ofile.printField("NN",nn); ofile.printField("MM",mm); 
  }
}

double NonLinearFunction::getCutoff( const double& width ) const {
  plumed_massert(setup,"non linear function has not been set");

  if(type==rational){
     plumed_merror("do not know what the cutoff is for rational functions");
  } else if(type==exponential){
     plumed_merror("do not know what the cutoff is for expoential functions");
  }  
 
  double DP2CUTOFF=6.25;
  if( type==gaussian ) return sqrt(2.0*DP2CUTOFF)*width;
  else if(type==triangular ) return width;
  else if(type==step) return width;
}

double NonLinearFunction::getVolume( const unsigned& dim, const double& determinant ) const {
  plumed_massert(setup,"non linear function has not been set");

  if(type==rational){
     plumed_merror("do not know what the volume is for rational functions");
  } else if(type==exponential){
     plumed_merror("do not know what the volume is for expoential functions");
  } 

  if( type==gaussian ){
     return pow( 2*pi, 0.5*dim ) * pow( determinant, 0.5 );
  } else {
     double vol;
     if( dim%2==1 ){
        double dfact=1;
        for(unsigned i=1;i<dim;i+=2) dfact*=static_cast<double>(i);
        vol=( pow( pi, (dim-1)/2 ) ) * ( pow( 2., (dim+1)/2 ) ) / dfact;
     } else {
        double fact=1.;
        for(unsigned i=1;i<dim/2;++i) fact*=static_cast<double>(i); 
        vol=pow( pi,dim/2 ) / fact;
     } 

     if(type==step) return vol*determinant;
     else if(type==triangular) return ( vol*determinant )/ 3.;
  }

  return 0;
}

