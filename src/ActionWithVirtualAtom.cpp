#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include "PlumedMain.h"

using namespace std;

namespace PLMD{

void ActionWithVirtualAtom::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.add("atoms","ATOMS","the list of atoms which are involved the virtual atom's definition");
}

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao)
{
  index=plumed.getAtoms().addVirtualAtom(this);
  AtomNumber a=AtomNumber::index(index);
  log.printf("  serial associated to this virtual atom is %d\n",a.serial());
}

ActionWithVirtualAtom::~ActionWithVirtualAtom(){
  plumed.getAtoms().removeVirtualAtom(this);
}

void ActionWithVirtualAtom::apply(){
  const Vector & f(plumed.getAtoms().forces[index]);
  for(unsigned i=0;i<getNumberOfAtoms();i++) modifyForces()[i]=matmul(derivatives[i],f);
}

void ActionWithVirtualAtom::requestAtoms(const std::vector<AtomNumber> & a){
  ActionAtomistic::requestAtoms(a);
  derivatives.resize(a.size());
}

void ActionWithVirtualAtom::setGradients(){
  Atoms&atoms(plumed.getAtoms());
  gradients.clear();
  for(unsigned i=0;i<getNumberOfAtoms();i++){
    AtomNumber an=getAbsoluteIndex(i);
    if(atoms.isVirtualAtom(an.index())){
      const ActionWithVirtualAtom* a=atoms.getVirtualAtomsAction(an.index());
      for(std::map<AtomNumber,Tensor>::const_iterator p=a->gradients.begin();p!=a->gradients.end();++p){
// controllare l'ordine del matmul:
        gradients[(*p).first]+=matmul(derivatives[i],(*p).second);
      }
    } else {
      gradients[an]+=derivatives[i];
    }
  }
}

}
