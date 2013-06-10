/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#ifndef __PLUMED_tools_Grid2_h
#define __PLUMED_tools_Grid2_h

#include <vector>
#include <string>
#include <cstring>
#include <iomanip>
#include <map>
#include <math.h>
#include "WeightBase.h"
// this below is to access to some testfiles without messing up the code
#include "config/Config.h"

namespace PLMD{ 


class Value;
class IFile;
class OFile;
class KernelFunctions;

/// \ingroup TOOLBOX
class Grid2
{
	std::vector<double> grid_;
	std::vector< std::vector<double> > der_;
protected:
	std::string funcname;
	std::vector<std::string> argnames;
	std::vector<std::string> str_min_, str_max_;
	std::vector<std::string> str_pmin_, str_pmax_;
	std::vector<double> min_,max_,dx_;
	std::vector<double> pmin_,pmax_,pdelta_,pshift_; // note that the pbc can be different from the min/max value
	std::vector<unsigned> nbin_,pbin_;
	std::vector<bool> pbc_;
	unsigned maxsize_, dimension_;
	bool dospline_, usederiv_;
	std::string fmt_; // format for output
	/// get "neighbors" for spline
	std::vector<unsigned> getSplineNeighbors(const std::vector<int> & indices)const;


public:
	/// clear grid
	virtual void clear();
	/// this constructor here is Value-aware
	///
	Grid2(const std::string& funcl, std::vector<Value*> args, const std::vector<std::string> & gmin,
			const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline,
			bool usederiv, bool doclear=true, bool docommensurate=true);
	/// this constructor here is not Value-aware
	Grid2(const std::string& funcl, const std::vector<std::string> &names, const std::vector<std::string> & gmin,
			const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline,
			bool usederiv, bool doclear, const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin,
			const std::vector<std::string> &pmax , bool docommensurate=true);
	/// this is the real initializator
	void Init(const std::string & funcl, const std::vector<std::string> &names, const std::vector<std::string> & gmin,
			const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline, bool usederiv,
			bool doclear, const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin, const std::vector<std::string> &pmax, bool docommensurate);
	/// get lower boundary
	std::vector<std::string> getMin() const;
	/// get upper boundary
	std::vector<std::string> getMax() const;
	/// get lower periodic boundary
	std::vector<std::string> getPeriodicMin() const;
	/// get upper periodic boundary
	std::vector<std::string> getPeriodicMax() const;
	/// get bin size
	std::vector<double> getDx() const;
	/// get bin volume
	double getBinVolume() const;
	/// get number of bins
	std::vector<unsigned> getNbin() const;
	/// get if periodic
	std::vector<bool> getIsPeriodic() const;
	/// get grid dimension
	unsigned getDimension() const;
	/// get argument names  of this grid
	std::vector<std::string> getArgNames() const;

	/// methods to handle grid indices
	/// indices can be positive or negative and are within the cell if the variable is periodic.
	/// return false when the index is not associated to any point
	/// return true if the conversion was successful
	bool getIndices(const unsigned index,std::vector<int> &) const;
	/// from value to a list of indices in each dimension
	/// can be positive or negative, depending on the location respect to min_
	/// if periodic: the index is always in the cell
	/// the conversion is always possible, therefore there is no return value
	void getIndices(const std::vector<double> & x, std::vector<int> &) const;
	// /// from a list of indices to a single index of the grid. If return bool is false then is out of grid
	bool getIndex(const std::vector<int> & indices, unsigned &index) const;
	// /// from a value to single index of the grid. if return bool is false is out of the grid
	bool getIndex(const std::vector<double> & x, unsigned &index) const;
	/// from a vector of values in full dimensionality to a vector of values right on the grid
	/// note that the conversion can always be done
	void getPoint(const std::vector<double> & x,std::vector<double> &xfloor) const;
	/// from an index to a point. Return bool is false whenever the index is out of grid
	bool getPoint(const unsigned index,std::vector<double> & point) const;
	/// from indices to a point in cv space. Conversion is always possible.
	void getPoint(const std::vector<int> & indices,std::vector<double> & point) const;
	/// when applied to a set of indices it returns the indices within the box, for periodic variables
	/// for nonperiodic variables the indices are untouched
	/// the return value is true if the value is within the grid and false elsewhere
	virtual  void applyPeriodicityToIndices(std::vector<int> &indices) const;
	// /// get neighbors
	std::vector<unsigned> getNeighbors(unsigned index,const std::vector<unsigned> & neigh) const;
	std::vector<unsigned> getNeighbors(const std::vector<int> & indices,const std::vector<unsigned> & neigh) const;
	std::vector<unsigned> getNeighbors(const std::vector<double> & x,const std::vector<unsigned> & neigh) const;
	//
	/// write header for grid file
	void writeHeader(OFile& file);

	/// static function that
	/// read grid from file
	/// and gives back the grid
	/// functl: the argumemnt to extract from the grid
	/// args: the arguments that are supposed to be there in the grid
	/// ifile: the input file that will be used to get the grid from
	/// dosparse: if you want to use sparse grids
	/// dospline: if you want to use spline
	/// doder: if you want to add derivative
	static Grid2* create(const std::string&,std::vector<Value*>,IFile&,bool dosparse,bool dospline,bool doder);
	// /// read grid from file and check boundaries are what is expected from input
	 static Grid2* create(const std::string&, std::vector<Value*>, IFile&,
			 const std::vector<std::string>&,const std::vector<std::string>&,
			 const std::vector<unsigned>&,bool,bool,bool);
	/// get grid size
	virtual unsigned getSize() const;
	// /// get grid value
	virtual bool getValue(unsigned index, double &val) const;
	virtual bool getValue(const std::vector<int> & indices, double &val) const;
	virtual bool getValue(const std::vector<double> & x, double &val ) const;
	/// get minimum value
	virtual double getMinValue() const;
	/// get maximum value
	virtual double getMaxValue() const;
	/// get grid value and derivatives
	virtual bool getValueAndDerivatives(unsigned index,double &val, std::vector<double>& der) const ;
	// gives false if you ask for a out-of-bound point of the grid
	virtual bool getValueAndDerivatives(const std::vector<int> & indices, double &val, std::vector<double>& der) const;
	// gives false if you ask for a out-of-bound point of the grid
	virtual bool getValueAndDerivatives(const std::vector<double> & x, double &val, std::vector<double>& der) const;
	//
	/// set grid value
	virtual bool setValue(const unsigned index, const double value);
	/// provide indices : return false if out of the grid
	virtual bool setValue(const std::vector<int> & indices, double value);
	/// set grid value and derivatives
	virtual bool setValueAndDerivatives(unsigned index, double value, std::vector<double>& der);
	/// return false when out of the grid
	virtual bool setValueAndDerivatives(const std::vector<int> & indices, double value, std::vector<double>& der);
	///return false when out of the grid
	virtual bool setValueAndDerivatives(const std::vector<double> & x, double value, std::vector<double>& der);
	/// add to grid value. return false if out of boundary
	virtual bool addValue(unsigned index, double value);
	//
	virtual bool addValue(const std::vector<int> & indices, double value);
	// /// add to grid value and derivatives
	virtual bool addValueAndDerivatives(unsigned index, double value, std::vector<double>& der);
	//
	virtual bool addValueAndDerivatives(const std::vector<int> & indices, double value, std::vector<double>& der);
	/// Scale all grid values and derivatives by a constant factor
	virtual void scaleAllValuesAndDerivatives( const double& scalef );
	// /// apply function: takes  pointer to  function that accepts a double and apply
	virtual void applyFunctionAllValuesAndDerivatives( double (*func)(double val), double (*funcder)(double valder) );
	/// add a kernel function to the grid
	void addKernel( const KernelFunctions& kernel );


	/// dump grid on file
	virtual void writeToFile(OFile&);

	virtual ~Grid2(){}

	/// project a high dimensional grid onto a low dimensional one: this should be changed at some time
	/// to enable many types of weighting
	/// proj: the values on which one should project
	/// ptr2obj: is a pointer to a class that allows to choose the suitable reweighting scheme
	/// gives back a new grid
	Grid2 project( const std::vector<std::string> & proj , WeightBase *ptr2obj  );
	void projectOnLowDimension(double &val , std::vector<int> &varHigh, WeightBase* ptr2obj );
	/// set output format
	void setOutputFmt(std::string ss){fmt_=ss;}
	/// gridtester. Try to test everything here
	static int gridTester();
};


class SparseGrid2 : public Grid2
{

 std::map<unsigned,double> map_;
 typedef std::map<unsigned,double>::const_iterator iterator;
 std::map< unsigned,std::vector<double> > der_;
 typedef std::map<unsigned,std::vector<double> >::const_iterator iterator_der;

 protected:
 void clear();

 public:
 SparseGrid2(const std::string& funcl, std::vector<Value*> args, const std::vector<std::string> & gmin,
            const std::vector<std::string> & gmax,
            const std::vector<unsigned> & nbin, bool dospline, bool usederiv,bool docommensurate=true):
            Grid2(funcl,args,gmin,gmax,nbin,dospline,usederiv,false,docommensurate){}

 // non value-aware
 SparseGrid2(const std::string& funcl, const std::vector<std::string> &names, const std::vector<std::string> & gmin,
            const std::vector<std::string> & gmax,
            const std::vector<unsigned> & nbin, bool dospline, bool usederiv, const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin,
			const std::vector<std::string> &pmax , bool docommensurate=true):
            Grid2(funcl,names,gmin,gmax,nbin,dospline,usederiv,false,isperiodic, pmin, pmax, docommensurate){}

 unsigned getSize() const;
 unsigned getMaxSize() const;

/// this is to access to Grid2:: version of these methods (allowing overloading of virtual methods)
 using Grid2::getValue;
 using Grid2::getValueAndDerivatives;
 using Grid2::setValue;
 using Grid2::setValueAndDerivatives;
 using Grid2::addValue;
 using Grid2::addValueAndDerivatives;

 /// get grid value
 bool getValue(unsigned index,double &val) const;
/// get grid value and derivatives
 bool getValueAndDerivatives(unsigned index, double &val, std::vector<double>& der) const;
/// set grid value
 bool setValue(unsigned index, double &value);
/// set grid value and derivatives
 bool setValueAndDerivatives(unsigned index, double &value, std::vector<double>& der);
/// add to grid value
 bool addValue(unsigned index, double &value);
/// add to grid value and derivatives
 bool addValueAndDerivatives(unsigned index, double value, std::vector<double>& der);
//
///// dump grid on file
 void writeToFile(OFile&);

 virtual ~SparseGrid2(){}
 };
}

#endif
