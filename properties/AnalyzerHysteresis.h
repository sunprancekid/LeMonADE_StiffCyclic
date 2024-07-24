/*--------------------------------------------------------------------------------
    ooo      L   attice-based  | LeMonADE: An Open Source Implementation of the
  o\.|./o    e   xtensible     |           Bond-Fluctuation-Model for Polymers
 o\.\|/./o   Mon te-Carlo      |           
oo---0---oo  A   lgorithm and  | ELMA-OscillatoryForce: Force-Extension-simulations
 o/./|\.\o   D   evelopment    | Copyright (C) 2022 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
    ooo      ELMA -            | Ron Dockhorn
             OscillatoryForce  |
----------------------------------------------------------------------------------

This file is part of LeMonADE and ELMA-OscillatoryForce extension.

LeMonADE and ELMA-OscillatoryForce extension is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE  and ELMA-OscillatoryForce extension is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE and ELMA-OscillatoryForce extension. If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef ANALYZERHYSTERESIS_H
#define ANALYZERHYSTERESIS_H

#include <vector>
#include <string>
#include <utility>      // std::pair
#include <map>          // std::map
#include <vector>
#include <algorithm>    // std::find_if
#include <iostream>
#include <functional>
#include <queue>
#include <cmath> //sqrt, abs
#include <fstream>
#include <sstream> // stringstream
#include <iomanip> // << setprecision (2) << fixed

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
//using namespace Eigen;


#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>



#include "StatisticMoment.h"

/* **************************************************************
 * A simple analyzer example for static properties 
 * 
 * In this case it calculates the average RgSquared of the molecules in
 * the system. Plain and simple...well...at least plain. Have fun.
 * *************************************************************/

template<class IngredientsType>
class AnalyzerHysteresis : public AbstractAnalyzer {
public:
    AnalyzerHysteresis(const IngredientsType& ing, std::string dstDir_, uint64_t evalulation_time_, uint64_t save_interval_);

    virtual ~AnalyzerHysteresis() {

    };

    //typedef typename IngredientsType::molecules_type molecules_type;
    const typename IngredientsType::molecules_type& molecules;

    const IngredientsType& getIngredients() const {
        return ingredients;
    }

    virtual void initialize();
    virtual bool execute();
    virtual void cleanup();
    
    bool isFileExisting(std::string);

private:

    const IngredientsType& ingredients;

    //only used to make sure you initialize your groups before you do things
    bool initialized;

    uint32_t counterFrames;

    std::string filename;
    std::string dstdir;

    // RNG
    RandomNumberGenerators rng;
    
    StatisticMoment Statistic_Rg2;

    StatisticMoment Statistic_Rx; // projected end-to-end distance only in x-direction (applied force)
    StatisticMoment Statistic_Ry; // projected end-to-end distance only in y-direction (not in force-direction)
    StatisticMoment Statistic_Rz; // projected end-to-end distance only in z-direction (not in force-direction)
    
    StatisticMoment Statistic_Ree; // end-to-end distance (Ree*Ree)^(1/2) only in  all direction
    
    StatisticMoment Statistic_Hysteresis; // average hysteresis area w/o absolute value
    StatisticMoment Statistic_HysteresisAbsolute; // average hysteresis area with absolute value
    
    //! Time in MCS, which need to pass to calculate the observables
    uint64_t evalulation_time;
    
    //! Time in MCS within one period, there the observables is collected
    uint64_t step_interval;
    
    uint32_t indexMonomersForce[2];
    
    //stores the projected end-to-end distance only in x-direction every savetime for one period T
    std::vector<double> AllRxxWithinPeriod;
    
    //stores the hysteresis integral after every period
    std::vector<double> vectorHysteresisIntegral;
    
    
    // calculation of hysteresis within on period using the trapezoidal method
    double calculateHysteresisInOnePeriodWithTrapezoidalMethod();
    
    // mapping t%T to the average projected end-to-end distance only in x-direction (applied force)
    std::map<uint64_t, StatisticMoment > timeToRx; 
    
    std::map<uint64_t, StatisticMoment >::iterator it;
};




/*****************************************************************************
 * IMPLEMENTATION OF METHODS
 * **************************************************************************/


/////////////////////////////////////////////////////////////////////////////

template<class IngredientsType>
AnalyzerHysteresis<IngredientsType>::AnalyzerHysteresis(const IngredientsType& ing, std::string dstDir_, uint64_t evalulation_time_, uint64_t save_interval_)
: ingredients(ing), molecules(ing.getMolecules()), dstdir(dstDir_), evalulation_time(evalulation_time_), step_interval(save_interval_) {
    counterFrames = 1;
    
    Statistic_Rg2.clear();
    
    Statistic_Rx.clear();
    Statistic_Ry.clear();
    Statistic_Rz.clear();
    
    Statistic_Ree.clear();
    
    Statistic_Hysteresis.clear();
    Statistic_HysteresisAbsolute.clear();
    
    AllRxxWithinPeriod.clear();
    
    vectorHysteresisIntegral.clear();
    
    initialized=false;
}

template<class IngredientsType>
void AnalyzerHysteresis<IngredientsType>::initialize() {
    
    //get the monomers with attributes subjected to force
	uint32_t n = 0;
	for (uint32_t i=0; i < ingredients.getMolecules().size();i++){
		if (ingredients.getMolecules()[i].getAttributeTag() != 0){
			indexMonomersForce[n] = i;
			n += 1;
		}
    }
		if((n > 2) || (n==0))
        {
            std::stringstream errormessage;
            errormessage<<"AnalyzerHysteresis::initialize()...ambigious number of monomers subjected to force" << std::endl; 
                    
            throw std::runtime_error(errormessage.str());
        }

    initialized=true;
}

template<class IngredientsType>
bool AnalyzerHysteresis<IngredientsType>::execute() {

    //check if groups have been initialized. if not, exit and explain
	if(initialized==false)
	{
		std::stringstream errormessage;
		errormessage<<"AnalyzerHysteresis::execute()...groups not initialized\n"
			    <<"Use AnalyzerHysteresis::initialize() or Taskmanager::init()\n";
			    
		throw std::runtime_error(errormessage.str());
	}
	
    if(ingredients.getMolecules().getAge() >= evalulation_time)
	{
        
        //Length of connection vector between the two charged Monomers = EndtoEnd
        double R_X = (ingredients.getMolecules()[indexMonomersForce[1]]).getX() - (ingredients.getMolecules()[indexMonomersForce[0]]).getX();
        double R_Y = (ingredients.getMolecules()[indexMonomersForce[1]]).getY() - (ingredients.getMolecules()[indexMonomersForce[0]]).getY();
        double R_Z = (ingredients.getMolecules()[indexMonomersForce[1]]).getZ() - (ingredients.getMolecules()[indexMonomersForce[0]]).getZ();
        
        VectorInt3 v3Ree = ingredients.getMolecules()[indexMonomersForce[1]] - ingredients.getMolecules()[indexMonomersForce[0]];
        
        Statistic_Rx.AddValue(R_X); // projected end-to-end distance only in x-direction (applied force)
        Statistic_Ry.AddValue(R_Y); // projected end-to-end distance only in y-direction (not in force-direction)
        Statistic_Rz.AddValue(R_Z); // projected end-to-end distance only in z-direction (not in force-direction)
    
        Statistic_Ree.AddValue(v3Ree.getLength()); // end-to-end distance (Ree*Ree)^(1/2) only in  all direction
    
        // order of instruction matters here!
        
        // calculate hysteresis int R(f) df for one period
        // check for beginning of new period and non-empty vector (first period) and delete all end-to-end vector information
        if ( ((ingredients.getMolecules().getAge()) % ingredients.getPeriodOscillatoryForce() == 0) && !AllRxxWithinPeriod.empty() ) {
            
            double hysteresis_integral = calculateHysteresisInOnePeriodWithTrapezoidalMethod();
            Statistic_Hysteresis.AddValue(hysteresis_integral);
            Statistic_HysteresisAbsolute.AddValue(std::abs(hysteresis_integral));
            
            vectorHysteresisIntegral.push_back(hysteresis_integral);
            
            AllRxxWithinPeriod.clear();
            
            std::cout << "run period: " << counterFrames << std::endl;
            counterFrames++;
        }
        
        // add recent projected end-to-end distance in x-direction
        AllRxxWithinPeriod.push_back(R_X);
        
        uint64_t modTime = (ingredients.getMolecules().getAge()) % ingredients.getPeriodOscillatoryForce();
        it = timeToRx.find(modTime);
        
        // already exist add the value
        if (it != timeToRx.end())
        {
            (it->second).AddValue(R_X);
        }
        else // create new value
        {
            timeToRx.insert(std::make_pair(modTime,StatisticMoment()));
            timeToRx[modTime].AddValue(R_X);
        }
        
    }
    
    return true;
}

// Using the trapezoidal method to integrate
// T_N (R) = 0.5*h*R(a) + 0.5*h*R(b) + sum_k=1^N-1 R(x_k)

template<class IngredientsType>
double AnalyzerHysteresis<IngredientsType>::calculateHysteresisInOnePeriodWithTrapezoidalMethod()
{
 double stepSize = step_interval;

 double omega = 2.0 * 3.14159265358979323846 / double(ingredients.getPeriodOscillatoryForce());
 double fA = ingredients.getAmplitudeOscillatoryForce();
 /* Finding Integration Value */
 
 // f = f0+fA*sin(omega*t) in eX-direction
 
 //         /         /T
 //        |          |
 // A(T) = o R * df = o R(f)*fA*omega*cos(omega*t) dt
 //        |          | 
 //        /         0/
 
 double hysteresis_integral = 0.0;
 
 for(size_t i = 1; i < AllRxxWithinPeriod.size(); i++)
 {
  hysteresis_integral += AllRxxWithinPeriod.at(i) * omega * fA *std::cos(omega * i*stepSize) * stepSize;
 }
 
 // The first value counts twice (2*0.5=1)at this is Rxx at t%T=0 in the cycle 
 // cos(omega * 0 * stepSize) = cos(omega * T * stepSize) != 1
 hysteresis_integral += AllRxxWithinPeriod.at(0) * omega * fA * stepSize;
 
 return hysteresis_integral;
 
}

/*
template<class IngredientsType>
double AnalyzerHysteresis<IngredientsType>::calculateHysteresisInOnePeriodWithSimpsonRule()
{
    
    stepSize = (upper - lower)/subInterval;
    
    
    integration = f(lower) + f(upper);
    
    for(i=1; i<= subInterval-1; i++)
    {
        k = lower + i*stepSize;
        
        if(i%2==0)
        {
            integration = integration + 2 * (f(k));
        }
        else
        {
            integration = integration + 4 * (f(k));
        }
        
    }
    
    integration = integration * stepSize/3;
}
*/
    
struct PathSeparator {

    bool operator()(char ch) const {
        return ch == '\\' || ch == '/';
    }
};

template<class IngredientsType>
bool AnalyzerHysteresis<IngredientsType>::isFileExisting(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

template<class IngredientsType>
void AnalyzerHysteresis<IngredientsType>::cleanup() {

    // print results into a file
    // get the filename and path
    // find the filename without path and extensions
    std::string filenameGeneral = std::string(std::find_if(ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator()).base(), ingredients.getName().end());

    std::string::size_type const p(filenameGeneral.find_last_of('.'));
    filenameGeneral = filenameGeneral.substr(0, p);

    std::stringstream systeminfo;
	systeminfo << "f0_" << std::setprecision (2) << std::fixed << ingredients.getOffsetOscillatoryForce() << "_"
               << "fA_" << std::setprecision (2) << std::fixed << ingredients.getAmplitudeOscillatoryForce() << "_"
               << "T_"  << ingredients.getPeriodOscillatoryForce();

    // construct a list
    std::vector < std::vector<double> > tmpResults_Hysteresis;

    // we have 3 columns and row
    uint32_t columns = 25;
    uint32_t rows = 1;



    // we have columns
    tmpResults_Hysteresis.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResults_Hysteresis[i].resize(rows);

   
    tmpResults_Hysteresis[0][0] = ingredients.getPeriodOscillatoryForce();
    tmpResults_Hysteresis[1][0] = 2.0*3.14159265358979323846/double(ingredients.getPeriodOscillatoryForce());Statistic_Rg2.ReturnM1();
    tmpResults_Hysteresis[2][0] = Statistic_Hysteresis.ReturnM1();
    tmpResults_Hysteresis[3][0] = Statistic_Hysteresis.ReturnM2();
    tmpResults_Hysteresis[4][0] = Statistic_Hysteresis.ReturnVar();
    tmpResults_Hysteresis[5][0] = 2.0 * Statistic_Hysteresis.ReturnSigma() / std::sqrt(1.0 * Statistic_Hysteresis.ReturnN());
    tmpResults_Hysteresis[6][0] = Statistic_Hysteresis.ReturnN();
    tmpResults_Hysteresis[7][0] = ingredients.getOffsetOscillatoryForce();
    tmpResults_Hysteresis[8][0] = ingredients.getAmplitudeOscillatoryForce();
    
    tmpResults_Hysteresis[9][0]  = Statistic_HysteresisAbsolute.ReturnM1();
    tmpResults_Hysteresis[10][0] = Statistic_HysteresisAbsolute.ReturnM2();
    tmpResults_Hysteresis[11][0] = Statistic_HysteresisAbsolute.ReturnVar();
    tmpResults_Hysteresis[12][0] = 2.0 * Statistic_HysteresisAbsolute.ReturnSigma() / std::sqrt(1.0 * Statistic_HysteresisAbsolute.ReturnN());
    tmpResults_Hysteresis[13][0] = Statistic_HysteresisAbsolute.ReturnN();
    
    
    tmpResults_Hysteresis[14][0] = Statistic_Ree.ReturnM1();
    tmpResults_Hysteresis[15][0] = Statistic_Ree.ReturnM2();
    tmpResults_Hysteresis[16][0] = Statistic_Ree.ReturnVar();
    tmpResults_Hysteresis[17][0] = 2.0 * Statistic_Ree.ReturnSigma() / std::sqrt(1.0 * Statistic_Ree.ReturnN());
    tmpResults_Hysteresis[18][0] = Statistic_Ree.ReturnN();

    tmpResults_Hysteresis[19][0] = Statistic_Rx.ReturnM1();
    tmpResults_Hysteresis[20][0] = Statistic_Rx.ReturnM2();
    tmpResults_Hysteresis[21][0] = Statistic_Ry.ReturnM1();
    tmpResults_Hysteresis[22][0] = Statistic_Ry.ReturnM2();
    tmpResults_Hysteresis[23][0] = Statistic_Rz.ReturnM1();
    tmpResults_Hysteresis[24][0] = Statistic_Rz.ReturnM2();
    

    std::stringstream comment;
    comment << " File produced by analyzer AnalyzerHysteresis" << std::endl
            << " for oscillatory force in x-direction: f = f0 + fA * sin(omega*t) " << std::endl
            << " calculating hysteresis A(omega) = circ int R(f) df = int_0^T R(f) * fA * omega * cos(omega*t) * dt" << std::endl
            << std::endl
            << std::endl
            << " T     ... period " << std::endl
            << " omega ... frequency = 2*PI/T " << std::endl
            << " f0    ... offset force " << std::endl
            << " fA    ... amplitude force " << std::endl
            << " <A>   ... hysteresis = <circ int R(f) df> " << std::endl
            << " <R>   ... end-to-end-distance = <(R*R)^(0.5)>= <(Rx^2+Ry^2+Rz^2)^(0.5)> " << std::endl
            << " <Rx>  ... end-to-end-vector component in x-direction" << std::endl
            << " <Ry>  ... end-to-end-vector component in y-direction" << std::endl
            << " <Rz>  ... end-to-end-vector component in z-direction" << std::endl
            << " under N = " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << std::endl
            << "T omega <A> <A²> Varianz(A) Error(A) SampleSize(A) f0 fA  <|A|> <|A|²> Varianz(|A|) Error(|A|) SampleSize(A) <R> <R²> Varianz(R) Error(R) SampleSize(R) <Rx> <Rx²> <Ry> <Ry²> <Rz> <Rz²> ";
            


    //new filename
    std::string filename_Hysteresis = filenameGeneral + "_" + systeminfo.str() + "_Hysteresis.dat";
    
    if(!isFileExisting(dstdir+"/"+filename_Hysteresis))
        ResultFormattingTools::writeResultFile(dstdir+"/"+filename_Hysteresis, this->ingredients, tmpResults_Hysteresis, comment.str());
    else
        ResultFormattingTools::appendToResultFile(dstdir + "/" + filename_Hysteresis, tmpResults_Hysteresis);

    
    
    // Average end-to-end-vector component Rx at times in period
    std::vector < std::vector<double> > tmpResults_Period;

    // we have columns and row
    columns = 7;
    rows = timeToRx.size();

    // we have columns
    tmpResults_Period.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResults_Period[i].resize(rows);
    
    // C++11 alternative:
    size_t counter = 0;
    
    for (const auto& entry : timeToRx) {
       
        double time = double(entry.first);
        
        double force = ingredients.getOffsetOscillatoryForce() + ingredients.getAmplitudeOscillatoryForce()*std::sin( (2.0*3.14159265359/     double(ingredients.getPeriodOscillatoryForce()) )*time);
        
        tmpResults_Period[0][counter] = time;
        tmpResults_Period[1][counter] = force;
        
        tmpResults_Period[2][counter] = entry.second.ReturnM1();
        tmpResults_Period[3][counter] = entry.second.ReturnM2();
        tmpResults_Period[4][counter] = entry.second.ReturnVar();
        tmpResults_Period[5][counter] = 2.0 * entry.second.ReturnSigma() / std::sqrt(1.0 * entry.second.ReturnN());
        tmpResults_Period[6][counter] = entry.second.ReturnN();
        
        counter++;
    }

    std::stringstream comment_Period;
    comment_Period << " File produced by analyzer AnalyzerHysteresis" << std::endl
            << " for oscillatory force in x-direction: f = f0 + fA * sin(omega*t) " << std::endl
            << " calculating hysteresis A(omega) = circ int R(f) df = int_0^T R(f) * fA * omega * cos(omega*t) * dt" << std::endl
            << std::endl
            << " T        ... period = " << ingredients.getPeriodOscillatoryForce() << std::endl
            << " omega    ... frequency = 2*PI/T = " << (2.0*3.14159265359/double(ingredients.getPeriodOscillatoryForce())) << std::endl
            << " f0       ... offset force = " << ingredients.getOffsetOscillatoryForce() << std::endl
            << " fA       ... amplitude force = " << ingredients.getAmplitudeOscillatoryForce()<< std::endl
            << " t        ... time " << std::endl
            << " <A>      ... hysteresis = <circ int R(f) df> " << std::endl
            << " <Rx>|t%T ... end-to-end-vector component in x-direction averaged at equal t%T over serveral periods" << std::endl
            << " under N = " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << std::endl
            << "t%T f(t%T) <Rx>|t%T <Rx²>|t%T Varianz(Rx) Error(Rx) SampleSize(Rx)";
           

    //new filename
    std::string filename_Period = filenameGeneral + "_" + systeminfo.str() + "_Period.dat";
    
    if(!isFileExisting(dstdir+"/"+filename_Period))
        ResultFormattingTools::writeResultFile(dstdir+"/"+filename_Period, this->ingredients, tmpResults_Period, comment_Period.str());
    else
        ResultFormattingTools::appendToResultFile(dstdir + "/" + filename_Period, tmpResults_Period);
    
    
    // ====================================================
    // List all entries of hysteresis area
    // ====================================================
    
    std::vector < std::vector<double> > tmpResults_HysteresisIntegralList;

    // we have columns and row
    columns = 1;
    rows = vectorHysteresisIntegral.size();

    // we have columns
    tmpResults_HysteresisIntegralList.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResults_HysteresisIntegralList[i].resize(rows);
    
    // C++11 alternative:
    size_t counter_vecHysteresisIntegral = 0;
    
    for (const auto& entry : vectorHysteresisIntegral) {
       
        tmpResults_HysteresisIntegralList[0][counter_vecHysteresisIntegral] = entry;
        
        counter_vecHysteresisIntegral++;
    }

    std::stringstream comment_HysteresisIntegral;
    comment_HysteresisIntegral << " File produced by analyzer AnalyzerHysteresis" << std::endl
            << " for oscillatory force in x-direction: f = f0 + fA * sin(omega*t) " << std::endl
            << " calculating hysteresis A(omega) = circ int R(f) df = int_0^T R(f) * fA * omega * cos(omega*t) * dt" << std::endl
            << " Prints every cicle of calculated hysteresis." << std::endl
            << std::endl
            << " T        ... period = " << ingredients.getPeriodOscillatoryForce() << std::endl
            << " omega    ... frequency = 2*PI/T = " << (2.0*3.14159265359/double(ingredients.getPeriodOscillatoryForce())) << std::endl
            << " f0       ... offset force = " << ingredients.getOffsetOscillatoryForce() << std::endl
            << " fA       ... amplitude force = " << ingredients.getAmplitudeOscillatoryForce()<< std::endl
            << " t        ... time " << std::endl
            << " A        ... hysteresis = circ int R(f) df " << std::endl
            << " under N = " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << std::endl
            << "A";
           

    //new filename
    std::string filename_HysteresisIntegral = filenameGeneral + "_" + systeminfo.str() + "_HysteresisIntegralList.dat";
    
    if(!isFileExisting(dstdir+"/"+filename_HysteresisIntegral))
        ResultFormattingTools::writeResultFile(dstdir+"/"+filename_HysteresisIntegral, this->ingredients, tmpResults_HysteresisIntegralList, comment_HysteresisIntegral.str());
    else
        ResultFormattingTools::appendToResultFile(dstdir + "/" + filename_HysteresisIntegral, tmpResults_HysteresisIntegralList);
    
    

}

#endif /*ANALYZERHYSTERESIS_H*/
