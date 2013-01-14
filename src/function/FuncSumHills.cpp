/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ActionRegister.h"
#include "Function.h"
#include "tools/Exception.h"
#include "tools/BiasRepresentation.h"
#include "tools/File.h"
#include "tools/Tools.h"
#include <iostream>

using namespace std;

namespace PLMD{
namespace function{


//+PLUMEDOC FUNCTION FUNCSUMHILLS 
/*

*/
//+ENDPLUMEDOC

class FilesHandler{
	vector <string> filenames; 
        vector <IFile*>  ifiles;
        Action *action;
        Log *log;
        bool parallelread;
        unsigned beingread;  
	public:
		FilesHandler(const vector<string> &filenames, const bool &parallelread ,  Action &myaction , Log &mylog);
		bool readBunch(BiasRepresentation *br, unsigned stride);
		bool scanOneHill(BiasRepresentation *br, IFile *ifile );
}; 
FilesHandler::FilesHandler(const vector<string> &filenames, const bool &parallelread , Action &action , Log &mylog ):filenames(filenames),parallelread(parallelread),beingread(0),log(&mylog){
   this->action=&action;
   for(unsigned i=0;i<filenames.size();i++){
      IFile *ifile = new IFile();
      ifile->link(action);
      ifiles.push_back(ifile);
      string ss;
      plumed_massert((ifile->FileExist(filenames[i])), "the file "+filenames[i]+" does not exist " );
   }
   
};

// note that the FileHandler is completely transparent respect to the biasrepresentation 
// no check are made at this level
bool FilesHandler::readBunch(BiasRepresentation *br , unsigned stride = -1){
	if(parallelread){
		(*log)<<"  doing parallelread \n";
        }else{
		(*log)<<"  doing serialread \n";
        	// read one by one hills      
		// is the type defined? if not, assume it is a gaussian 
                IFile *ff; 
                ff=ifiles[beingread];
                ff->open(filenames[beingread]);
		unsigned n=0;
                
                while(n<stride && stride >0 ){
			while(scanOneHill(br,ff)){
 				n++;
                        }   // read one hill and update the bias representation
                        (*log)<<"  closing file "<<filenames[beingread]<<"\n";
                        ff->close();
			beingread++;
 		        if(beingread<ifiles.size()){ff=ifiles[beingread];ff->open(filenames[beingread]);}else{break;}  
                } 
        }        
	return true;
};
bool FilesHandler::scanOneHill(BiasRepresentation *br, IFile *ifile ){
	// return false if the file is over	
	std::string kernel;
	double dummy;
	if(ifile->scanField("time",dummy)){
	        (*log)<<"   scanning one hill: "<<dummy<<" \n";
                // find the kernel
	        if(ifile->FieldExist("kernel")){ifile->scanField("kernel",kernel);
	        }else{
	          (*log)<<"  kernel is assumed Gaussian \n"; 
                  kernel="gaussian";
	        };	 
                // scan the line
                unsigned ncv; 
                ncv=br->getNumberOfDimensions();
                vector<double> center(ncv); 
                Value* v;
		for(unsigned i=0;i<ncv;i++){
			v=br->getPtrToValue(i);
 		 	(*log)<<"scanning "<<br->getName(i)<<"\n";	
			ifile->scanField(v);
			center[i]=v->get();		
			//ifile->scanField( br->getName(i), center[i]);
	        	std::string imin, imax;
                	v->getDomain( imin, imax );
		}
		// kernel dependent reading
	        if(kernel=="gaussian"){
			bool multivariate;
			std::string sss;
			ifile->scanField("multivariate",sss);
			if(sss=="true") multivariate=true;
			else if(sss=="false") multivariate=false;
			else plumed_merror("cannot parse multivariate = "+ sss);
                        vector<double> sigma; 
			if(multivariate){
			   sigma.resize(ncv*(ncv+1)/2);
			   Matrix<double> upper(ncv,ncv);
			   Matrix<double> lower(ncv,ncv);
			   for (unsigned i=0;i<ncv;i++){
			            for (unsigned j=0;j<ncv-i;j++){
                                            string mysigma;
					    mysigma="sigma_"+(br->getPtrToValue(j+i))->getName()+"_"+(br->getPtrToValue(j))->getName();
			                    //ifile->scanField("sigma_"+(br->getValue(j+i))->getName()+"_"+(br->getValue(j))->getName(),lower(j+i,j));
			                    ifile->scanField(mysigma,lower(j+i,j));
			                    upper(j,j+i)=lower(j+i,j);
			            }
			   }
			   Matrix<double> mymult(ncv,ncv);       
			   Matrix<double> invmatrix(ncv,ncv);       
			   //log<<"Lower \n";
			   //matrixOut(log,lower); 
			   //log<<"Upper \n";
			   //matrixOut(log,upper); 
			   mult(lower,upper,mymult);          
			   //log<<"Mult \n";
			   //matrixOut(log,mymult); 
			   // now invert and get the sigmas
			   Invert(mymult,invmatrix);
			   //log<<"Invert \n";
			   //matrixOut(log,invmatrix); 
			   // put the sigmas in the usual order: upper diagonal (this time in normal form and not in band form) 
			   unsigned k=0;
			   for (unsigned i=0;i<ncv;i++){
			          for (unsigned j=i;j<ncv;j++){
			          	sigma[k]=invmatrix(i,j);
			          	k++;
			          }
			   }
			}else{
				for(unsigned i=0;i<ncv;++i){
			       ifile->scanField("sigma_"+(br->getPtrToValue(i))->getName(),sigma[i]);
			   }
			}
			double height,dummy;
			ifile->scanField("height",height);
			ifile->scanField("biasf",dummy);
			if(ifile->FieldExist("clock")) ifile->scanField("clock",dummy);
                	// now push this into the bias representation
			// the centers are already inside the updated values
                        br->pushGaussian( sigma, height);
                }else{
			plumed_merror("This kernel reader is not implemented");
                };	
	 	(*log)<<"  read hill\n";	
	        ifile->scanField();  
 		return true;
	}else{
		return false;
	}
};

class FuncSumHills :
  public Function
{
  vector<string> hillsFiles,histoFiles; 
  vector<string> proj; 
  unsigned initstride,stride,highdim, lowdim;
  bool iscltool,integratehills,integratehisto,parallelread;
  double beta;
  BiasRepresentation *biasrep;
  BiasRepresentation *historep;
public:
  FuncSumHills(const ActionOptions&);
  ~FuncSumHills();
  void calculate(); // this probably is not needed
  bool checkFilesAreExisting(const vector<string> hills ); 
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FuncSumHills,"FUNCSUMHILLS")

void FuncSumHills::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); 
  keys.add("optional","HILLSFILES"," source file for hills creation(may be the same as HILLS)"); // this can be a vector! 
  keys.add("optional","HISTOFILES"," source file for histogram creation(may be the same as HILLS)"); // also this can be a vector!
  keys.add("optional","PROJ"," only with sumhills: the projection on the cvs");
  keys.add("optional","KT"," only with sumhills: the kt factor when projection on cvs");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid"); 
  keys.add("optional","OUTHILLS"," output file for hills ");
  keys.add("optional","INITSTRIDE"," stride if you want an initial dump ");
  keys.add("optional","STRIDE"," stride when you do it on the fly ");
  keys.addFlag("ISCLTOOL",true,"use via plumed commandline: calculate at read phase and then go");
  keys.addFlag("PARALLELREAD",false,"read parallel HILLS file");
}

FuncSumHills::FuncSumHills(const ActionOptions&ao):
Action(ao),
Function(ao),
initstride(-1),
stride(-1),
iscltool(false),
beta(-1.),
integratehills(true),
integratehisto(false),
parallelread(false)
{
  // here read 
  // Grid Stuff
  vector<std::string> gmin(getNumberOfArguments());
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments() && gmin.size()!=0) error("not enough values for GRID_MIN");
  plumed_massert(gmin.size()==getNumberOfArguments(),"need GRID_MIN argument for this") ;
  vector<std::string> gmax(getNumberOfArguments());
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments() && gmax.size()!=0) error("not enough values for GRID_MAX");
  plumed_massert(gmax.size()==getNumberOfArguments(),"need GRID_MAX argument for this") ;
  vector<unsigned> gbin(getNumberOfArguments());
  parseVector("GRID_BIN",gbin);
  plumed_massert(gbin.size()==getNumberOfArguments(),"need GRID_BIN argument for this"); 
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("not enough values for GRID_BIN");
  plumed_assert(gmin.size()==gmax.size() && gmin.size()==gbin.size());
  // add some automatic hills width: not in case stride is defined  
  // since when you start from zero the automatic size will be zero!
 
  // hills file: 
  parseVector("HILLSFILES",hillsFiles);
  if(hillsFiles.size()==0){
  	hillsFiles.push_back("HILLS");
  	integratehills=true; // default behaviour  
  }
  // histo file: 
  parseVector("HISTOFILES",histoFiles);
  if(histoFiles.size()==0){
  	integratehisto=false;  
  }
  // needs a projection? 
  proj.clear();
  parseVector("PROJ",proj);
  plumed_massert(proj.size()<getNumberOfArguments()," The number of projection must be less than the full list of arguments ");
  string kt="";
  if(proj.size()>0 ) {
    parse("KT",kt);
    plumed_massert(kt!="","if you make a projection then you need KT flag!"); 
    Tools::convert(beta,kt); beta=1./beta; 
  }
  // is a cltool: then you start and then die
  parseFlag("ISCLTOOL",iscltool);
  //
  parseFlag("PARALLELREAD",parallelread);
  // stride
  string initstride_s;
  parse("INITSTRIDE",initstride_s);
  Tools::convert(initstride,initstride_s);

  //what might it be this? 
  checkRead();
  // here start 
  // want something right now?? do it and return
  // your argument is a set of cvs 
  // then you need: a hills / a colvar-like file (to do a histogram) 
  // create a bias representation for this
  if(iscltool){

    //std::vector<Value> tmpvalues; 
    std::vector<Value*> tmpvalues; 
    for(unsigned i=0;i<getNumberOfArguments();i++){
        // allocate a new value from the old one: no deriv here
	//tmpvalues.push_back(  Value( this, getPntrToArgument(i)->getName(), true ) );
	tmpvalues.push_back( getPntrToArgument(i) );
    }

    // check if the files exists 
    if(integratehills){
         checkFilesAreExisting(hillsFiles); 
         biasrep=new BiasRepresentation(tmpvalues,comm, gmax,gmin,gbin);
    }
    if(integratehisto){
         checkFilesAreExisting(histoFiles); 
         historep=new BiasRepresentation(tmpvalues,comm,gmax,gmin,gbin);
    }

    // decide how to source hills ( serial/parallel )
    // here below the input control 
    // say how many hills and it will read them from the 
    // bunch of files provided, will update the representation 
    // of hills (i.e. a list of hills and the associated grid)

    // decide how to source colvars ( serial parallel )
    FilesHandler *hillsHandler=new FilesHandler(hillsFiles,parallelread,*this, log);

    // read a number of hills and put in the bias representation
    hillsHandler->readBunch(biasrep,initstride);


    // project and dump

    return;
  } 
  // just an initialization but you need to do something on the fly?: need to connect with a metad run and its grid representation 
  // your argument is a metad run
  // if the grid does not exist crash and say that you need some data 
  // otherwise just link with it

}

void FuncSumHills::calculate(){
  // this should be connected only with a grid representation to metadynamics 
  // at regular time just dump it
}

FuncSumHills::~FuncSumHills(){
}

bool FuncSumHills::checkFilesAreExisting(const vector<string> hills ){
	plumed_massert(hills.size()!=0,"the number of  files provided should be at least one" );
        IFile *ifile = new IFile();
        ifile->link(*this);
        for(unsigned i; i< hills.size();i++){  
          plumed_massert(ifile->FileExist(hills[i]),"missing file "+hills[i]);
        }
        return true; 

}

}

}


