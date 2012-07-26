#include "RMSDMap.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD {

void RMSDMap::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys ); 
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys ); 
  keys.add("numbered","REFERENCE","a file in pdb format containing structure at one of the reference point. " + PDB::documentation() ); 
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.add("compulsory","LAMBDA","the lamdba parameter for the path");
  keys.reserve("numbered","LOW_DIM_VEC","the projections to use for each of the reference points"); 
  ActionWithDistribution::registerKeywords( keys );
}

RMSDMap::RMSDMap(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithDistribution(ao)
{

 std::string type; type.assign("SIMPLE"); parse("TYPE",type);
 std::vector<AtomNumber> save_atoms; 
 unsigned nldim; std::vector<double> vecl;
 for(unsigned i=1;;++i){
    std::string reference; parseNumbered("REFERENCE",i,reference);
    if( reference.length()==0 ) break;
    // Read the pdb file 
    PDB pdb;
    if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().length) )
         error("missing input file " + reference);
    std::vector<AtomNumber> atoms(pdb.getAtomNumbers());
    if(i==1){
       requestAtoms(atoms); save_atoms.resize( atoms.size() );
       for(unsigned i=0;i<atoms.size();++i) save_atoms[i]=atoms[i];
    } else {
       if(atoms.size()!=save_atoms.size() ) error("mismatched atoms in pdb files");
       for(unsigned i=0;i<atoms.size();++i){
           if( atoms[i]!=save_atoms[i] ) error("mismatched atoms in pdb files");
       }
    }
    // Save the positions of the atoms
    frames.push_back( RMSD() ); 
    frames[i-1].set(pdb,type); 
    log.printf("  projecting configuration in %s ",reference.c_str() );
    // Read in the low dimensional vectors
    if( keywords.exists("LOW_DIM_VEC") ){
        parseNumberedVector("LOW_DIM_VEC",i,vecl);
        if( i!=1 && vecl.size()!=nldim ) error("size mismatches for projections");
        else if( i==1 ) nldim=vecl.size();
        log.printf("at : ");
        for(unsigned i=0;i<vecl.size();++i) log.printf("%f ",vecl[i]);
        log.printf("\n"); 
        low_dims.push_back( vecl );
    } else {
        double pp=static_cast<double>(i);
        log.printf("at : %f\n",pp);
        low_dims.push_back( std::vector<double>( 1, pp ) );
    } 
 }
 // Read the lambda parameter
 parse("LAMBDA",lambda);

 // Read distribution keywords and resize the functions
 requestDistribution(); 

 // Resize everything
 derivs.resize( getNumberOfAtoms() );
 thevalue.resizeDerivatives( 3*getNumberOfAtoms() + 9 );
}

void RMSDMap::calculate(){
  calculateAllVessels( getStep() );
}

bool RMSDMap::calculateThisFunction( const unsigned& j ){
  plumed_assert( j<frames.size() );
 
  thevalue.clearDerivatives(); vir.zero();
  double r=frames[j].calculate(getPositions(),derivs,log,true); 
  for(unsigned i=0;i<derivs.size();i++) vir=vir+(-1.0*Tensor(getPosition(i),derivs[i]));

  thevalue.set( exp(-r*lambda ) ); double der=-lambda*exp(-r*lambda );
  for(unsigned i=0;i<derivs.size();i++){
      thevalue.addDerivative( 3*i+0, der*derivs[i][0] );
      thevalue.addDerivative( 3*i+1, der*derivs[i][1] );
      thevalue.addDerivative( 3*i+2, der*derivs[i][2] ); 
  }
  unsigned natoms=getNumberOfAtoms();
  thevalue.addDerivative( 3*natoms + 0, der*vir(0,0) );
  thevalue.addDerivative( 3*natoms + 1, der*vir(0,1) );
  thevalue.addDerivative( 3*natoms + 2, der*vir(0,2) );
  thevalue.addDerivative( 3*natoms + 3, der*vir(1,0) );
  thevalue.addDerivative( 3*natoms + 4, der*vir(1,1) );
  thevalue.addDerivative( 3*natoms + 5, der*vir(1,2) );
  thevalue.addDerivative( 3*natoms + 6, der*vir(2,0) );
  thevalue.addDerivative( 3*natoms + 7, der*vir(2,1) );
  thevalue.addDerivative( 3*natoms + 8, der*vir(2,2) );
  return false;
}

void RMSDMap::apply(){
  std::vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());

  for(unsigned i=0;i<f.size();i++){
    f[i][0]=0.0;
    f[i][1]=0.0;
    f[i][2]=0.0;
  }
  v.zero();

  unsigned nat=getNumberOfAtoms();
  std::vector<double> forces(3*getNumberOfAtoms()+9);

  unsigned vstart=3*getNumberOfAtoms();
  for(int i=0;i<getNumberOfVessels();++i){
    if( (getPntrToVessel(i)->applyForce( forces )) ){
     for(unsigned j=0;j<nat;++j){
        f[j][0]+=forces[3*j+0];
        f[j][1]+=forces[3*j+1];
        f[j][2]+=forces[3*j+2];
     }
     v(0,0)+=forces[vstart+0];
     v(0,1)+=forces[vstart+1];
     v(0,2)+=forces[vstart+2];
     v(1,0)+=forces[vstart+3];
     v(1,1)+=forces[vstart+4];
     v(1,2)+=forces[vstart+5];
     v(2,0)+=forces[vstart+6];
     v(2,1)+=forces[vstart+7];
     v(2,2)+=forces[vstart+8];
    }
  }
}

}


