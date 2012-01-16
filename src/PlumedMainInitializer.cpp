#include "PlumedMain.h"
#include "PlumedException.h"
#include <cstdlib>

using namespace PLMD;
using namespace std;

// !!!!!!!!!!!!!!!!!!!!!!    DANGER   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// THE FOLLOWING ARE UTILITIES WHICH ARE NECESSARY FOR DYNAMIC LOADING OF THE PLUMED KERNEL:
// This section should be consistent with the Plumed.h file.
// Since the Plumed.h file may be included in host MD codes, **NEVER** MODIFY THE CODE IN THIS FILE
// unless you know exactly what you are doing

/* Holder for plumedmain function pointers */
typedef struct {
  void*(*create)();
  void(*cmd)(void*,const char*,const void*);
  void(*finalize)(void*);
} plumed_plumedmain_function_holder;

/* These functions should be accessible from C, since they might be statically
   used from Plumed.c (for static binding) */

extern "C" void*plumedmain_create(){
  return new PlumedMain;
}

extern "C" void plumedmain_cmd(void*plumed,const char*key,const void*val){
  plumed_massert(plumed,"trying to use a plumed object which is not initialized");
  static_cast<PlumedMain*>(plumed)->cmd(key,val);
}

extern "C" void plumedmain_finalize(void*plumed){
  plumed_massert(plumed,"trying to deallocate a plumed object which is not initialized");
  delete static_cast<PlumedMain*>(plumed);
}

/* This refers to a function implemented in Plumed.c */
extern "C" plumed_plumedmain_function_holder* plumed_kernel_register(const plumed_plumedmain_function_holder*);

namespace PLMD{

/// Static object which registers Plumed.
/// This is a static object which, during its construction at startup,
/// registers the pointers to plumedmain_create, plumedmain_cmd and plumedmain_finalize
/// to the plumed_kernel_register function
static class PlumedMainInitializer{
  public:
  PlumedMainInitializer(){
    plumed_plumedmain_function_holder fh={plumedmain_create,plumedmain_cmd,plumedmain_finalize};
    plumed_kernel_register(&fh);
  };
} RegisterMe;

}

