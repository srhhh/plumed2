#include "Value.h"
#include "ActionWithValue.h"
#include "PlumedException.h"

using namespace PLMD;

Value::Value(const std::string& name, const bool withderiv):
  value(0.0),
  name(name),
  hasDeriv(withderiv),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{
}

void Value::setupPeriodicity(){
  if( min==0 && max==0 ){
      periodicity=notperiodic;
  } else {
      periodicity=periodic;
      max_minus_min=max-min;
      plumed_massert(max_minus_min>0, "your function has a very strange domain?");
      inv_max_minus_min=1.0/max_minus_min;
  }
}

bool Value::isPeriodic()const{
  plumed_massert(periodicity!=unset,"periodicity should be set");
  return periodicity==periodic;
}

bool Value::applyForce(std::vector<double>& forces ) const {
  plumed_massert( derivatives.size()==forces.size()," forces array has wrong size" );
  if( !hasForce ) return false;
  for(unsigned i=0;i<derivatives.size();++i) forces[i]=inputForce*derivatives[i]; 
  return true;
}

void Value::getDomain(double&min,double&max) const {
  plumed_massert(periodicity==periodic,"function should be periodic");
  min=this->min;
  max=this->max;
}
