#include "FieldBias.h"
#include "PlumedMain.h"
#include "ActionSet.h"

using namespace std;
using namespace PLMD;

void FieldBias::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which to output the field");
  keys.add("compulsory","FIELD","the input for this action is the field calculated during one of the other actions.");
  keys.add("compulsory","NORM","the normalization to use");
  keys.add("optional","NGRID","number of grid points to use in each direction - if you are using function interpolation");
  keys.add("optional","START_BIAS","the bias at the start of the simulation");
  keys.addFlag("SERIAL", false, "do the calculation in serial");
  keys.addFlag("DEBUG_DERIVATIVES",false,"used to debug the derivatives of the bias");
}

FieldBias::FieldBias(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionPilot(ao),
  serial(false),
  debug(false)
{
  if( checkNumericalDerivatives() ) warning("cannot parallelize field with numerical derivatives over actions");
  parseFlag("SERIAL",serial);
  if(serial) log.printf("  doing calculation in serial\n");
  parseFlag("DEBUG_DERIVATIVES",debug);
  if(debug) log.printf("  storing derivatives in a value so they can be output with dump derivatives");

  // Find the field we are using
  std::string ll; parse("FIELD",ll);
  ActionWithDistribution* field=plumed.getActionSet().selectWithLabel<ActionWithDistribution*>(ll);
  addDependency(field);
  if(!field) error("cannot find action named " + ll);
  myfield=field->getField();
  if(!myfield) error("action " + ll + " calculates a colvar and not a field");
  apply_action=dynamic_cast<ActionWithDistribution*>( field );

  // Read how we descritize the integrals
  std::vector<unsigned> ngrid( myfield->get_Ndx() ); 
  parseVector("NGRID",ngrid);
  if( ngrid.size()==0 ){ ngrid.resize( myfield->get_Ndx() ); myfield->get_nspline( ngrid ); }
  if( ngrid.size()!=myfield->get_Ndx() ) error("grid size is wrong");
  log.printf("  using field %s and descritizing %d dimensional integrals over ",ll.c_str(), ngrid.size() );
  for(unsigned i=0;i<ngrid.size();++i) log.printf("%d ",ngrid[i]);
  // Work out the norm we are using
  unsigned nn; parse("NORM",nn); norm=static_cast<double>(nn);
  log.printf("points. Normalizing with %d norm\n",nn);

  // Create the grid where we store the bias
  std::vector<double> min, max;
  myfield->retrieveBoundaries( min, max );
  plumed_assert( min.size()==ngrid.size() && max.size()==ngrid.size() );

  std::string sbias; parse("START_BIAS",sbias);
  if( sbias.length()==0 ){
    std::vector<bool> pbc(min.size(), false );
    bias=new Grid( min, max, ngrid, pbc, false, false );
  } else {
    log.printf("  reading initial bias from file named %s\n",sbias.c_str());
    FILE* bfile=fopen( sbias.c_str(),"r");
    bias=Grid::create( bfile, false, false, false );
    fclose( bfile );
    for(unsigned i=0;i<ngrid.size();++i){
       if( bias->getMin()[i]!=min[i] ) error("minimum in input grid does not match minimum of field");
       if( ( bias->getMax()[i]- (max[i] + bias->getDx()[i]) )>0.00001 ) error("maximum in input grid does not match maximum of field");
       if( bias->getNbin()[i]!=(ngrid[i]+1) ) error("number of bins in input grid does not match plumed input");
       if( bias->getIsPeriodic()[i]==true ) error("periodic input grid is not compatible with field overlap");
    }
  }
  // Prepare the buffers for the calculation
  buffer.resize( bias->getSize() + 2 );

  // Setup the blocks for parallelizing the grid calculations
  unsigned stride=comm.Get_size();
  if(serial) stride=1;

  blocks.resize( stride+1 );
  nn=std::floor( bias->getSize() / stride );
  unsigned nrem=bias->getSize() - nn*stride;

  blocks[0]=0;
  for(unsigned i=1;i<blocks.size();++i){
      for(unsigned j=0;j<i;++j) blocks[i]+=blocks[j];
      if( i<=nrem ) blocks[i]+=nn + 1;
      else blocks[i]+=nn;
  }
  plumed_assert( blocks[blocks.size()-1]==bias->getSize() );

  // Add something for the bias and set up the forces
  addComponentWithDerivatives("bias"); 
  if( debug ) getPntrToComponent("bias")->resizeDerivatives( 1 );
}

void FieldBias::clearBias(){
  std::vector<unsigned> ngrid; ngrid=bias->getNbin();
  delete bias;
  std::vector<double> min, max; myfield->retrieveBoundaries( min, max );
  plumed_assert( min.size()==ngrid.size() && max.size()==ngrid.size() );
  std::vector<bool> pbc(min.size(), false );
  bias=new Grid( min, max, ngrid, pbc, false, false );
}

void FieldBias::calculate(){
  unsigned rank=comm.Get_rank();
  if(serial) rank=0;

  unsigned nder=myfield->get_NdX();
  if( derivatives.size()!=nder ){
      derivatives.resize( nder ); 
      if( debug ) getPntrToComponent("bias")->resizeDerivatives( apply_action->getNumberOfDerivatives() ); 
  }
  // This calculates the current field if we are doing numerical derivatives
  if( checkNumericalDerivatives() ) apply_action->calculate(); 

  // The loop for the bias
  buffer.assign( buffer.size(), 0.0 );
  std::vector<double> pp( myfield->get_Ndx() );
  for(unsigned i=blocks[rank];i<blocks[rank+1];++i){
      bias->getPoint( i, pp );
      buffer[i+2]=myfield->calculateField( pp );
      buffer[0]+=pow(buffer[i+2], norm);
      buffer[1]+=buffer[i+2]*bias->getValue( i );
  }
  buffer[0]*=bias->getBinVolume();   
  if(!serial) comm.Sum( &buffer[0], buffer.size() );
  double normali=pow( buffer[0], 1./static_cast<double>(norm) );
  buffer[1]*=bias->getBinVolume() / normali;
  getPntrToComponent("bias")->set( buffer[1]  );

  if( checkNumericalDerivatives() ) return ;
 
  // The loop for the derivatives 
  derivatives.assign( derivatives.size(), 0.0 );
  std::vector<double> tmpforce( derivatives.size() );
  for(unsigned i=blocks[rank];i<blocks[rank+1];++i){
      bias->getPoint( i, pp );
      myfield->calculateFieldDerivatives( pp, tmpforce ); 
      for(unsigned j=0;j<derivatives.size();++j){
          derivatives[j] += tmpforce[j] * ( bias->getValue( i ) / normali - (buffer[1]/buffer[0])*pow( buffer[i+2], norm-1) ); 
      }
  }
  for(unsigned j=0;j<derivatives.size();++j) derivatives[j]*=bias->getBinVolume();
  if(!serial) comm.Sum( &derivatives[0], derivatives.size() );
  if( debug ) apply_action->mergeFieldDerivatives( derivatives, getPntrToComponent("bias") );
}

void FieldBias::calculateNumericalDerivatives( ActionWithValue* a ){
  apply_action->calculateNumericalDerivatives( this );
}

void FieldBias::apply(){
  if(onStep()){
     for(unsigned j=0;j<derivatives.size();++j) derivatives[j]*=-1.0*getStride();
     myfield->addForces( derivatives ); 
  }
}

FieldBias::~FieldBias(){
  delete bias; 
}


