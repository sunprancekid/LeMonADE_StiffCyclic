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

#ifndef ANALYZER_CONSTANT_FORCE_H
#define ANALYZER_CONSTANT_FORCE_H

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

// packages not included in LeMonADE
#include "StatisticMoment.h"

/* **************************************************************
 * A simple analyzer example for static properties 
 * 
 * In this case it calculates the average RgSquared of the molecules in
 * the system. Plain and simple...well...at least plain. Have fun.
 * *************************************************************/

template<class IngredientsType>
class AnalyzerConstantForce : public AbstractAnalyzer {
public:
    AnalyzerConstantForce(const IngredientsType& ing, std::string dstDir_, uint64_t evalulation_time_);

    virtual ~AnalyzerConstantForce() {

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
    
    StatisticMoment Statistic_Rg2; // radius of gyration 
    StatisticMoment Statistic_Rg2_x; // radius of gyration in x
    StatisticMoment Statistic_Rg2_y; // radius of gyration in y
    StatisticMoment Statistic_Rg2_z; // radius of gyration in z

    StatisticMoment Statistic_Rx; // projected end-to-end distance only in x-direction (applied force)
    StatisticMoment Statistic_Ry; // projected end-to-end distance only in y-direction (not in force-direction)
    StatisticMoment Statistic_Rz; // projected end-to-end distance only in z-direction (not in force-direction)
    
    StatisticMoment Statistic_Ree; // end-to-end distance (Ree*Ree)^(1/2) only in  all direction
    
    //! Time in MCS, which need to pass to calculate the observables
    uint64_t evalulation_time;
    
    
    uint32_t indexMonomersForce[2];
    
    //stores the end-to-end-vector component in x-direction
    std::vector<double> vectorAllRxx;

};

/*****************************************************************************
 * IMPLEMENTATION OF METHODS
 * **************************************************************************/

/////////////////////////////////////////////////////////////////////////////

template<class IngredientsType>
AnalyzerConstantForce<IngredientsType>::AnalyzerConstantForce(const IngredientsType& ing, std::string dstDir_, uint64_t evalulation_time_)
: ingredients(ing), molecules(ing.getMolecules()), dstdir(dstDir_), evalulation_time(evalulation_time_) {
    counterFrames = 1;
    
    Statistic_Rg2.clear();
    Statistic_Rg2_x.clear();
    Statistic_Rg2_y.clear();
    Statistic_Rg2_z.clear();
    
    Statistic_Rx.clear();
    Statistic_Ry.clear();
    Statistic_Rz.clear();
    
    Statistic_Ree.clear();
   
    vectorAllRxx.clear();
    
    initialized=false;
}

template<class IngredientsType>
void AnalyzerConstantForce<IngredientsType>::initialize() {
    
    //get the monomers with attributes subjected to force
	uint32_t n = 0;
	for (uint32_t i=0; i < ingredients.getMolecules().size();i++){
		if (ingredients.getMolecules()[i].getAttributeTag() != 0){
            std::cout << ingredients.getMolecules()[i].getAttributeTag() << "\n";
			indexMonomersForce[n] = i;
			n += 1;
		}
    }

    if((n > 2) || (n==0))
    {
        std::stringstream errormessage;
        errormessage << n << "\n" << std::endl;
        errormessage<<"AnalyzerConstantForce::initialize()...ambigious number of monomers subjected to force" << std::endl;

        throw std::runtime_error(errormessage.str());
    }

    initialized=true;
}

template<class IngredientsType>
bool AnalyzerConstantForce<IngredientsType>::execute() {

    //check if groups have been initialized. if not, exit and explain
	if(initialized==false)
	{
		std::stringstream errormessage;
		errormessage<<"AnalyzerConstantForce::execute()...groups not initialized\n"
			    <<"Use AnalyzerConstantForce::initialize() or Taskmanager::init()\n";
			    
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
    
        vectorAllRxx.push_back(R_X); // projected end-to-end vector only in x-direction (applied force)
        
        // The radius of gyration is defined as
		// Rg2 = 1/N * SUM_i=1 (r_i - r_COM)^2 = 1/N^2 * SUM_i=1 SUM_j=i (r_i - r_j)^2
        int monomerCounter = 0;
		double Rg2 = 0.0;
		double Rg2_x = 0.0;
		double Rg2_y = 0.0;
		double Rg2_z = 0.0;

		for (int k= 0; k < ingredients.getMolecules().size(); k++)
		{
			for (int l= k; l < ingredients.getMolecules().size(); l++)
			{
				Rg2_x += (ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX())*(ingredients.getMolecules()[k].getX()-ingredients.getMolecules()[l].getX());
				Rg2_y += (ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY())*(ingredients.getMolecules()[k].getY()-ingredients.getMolecules()[l].getY());
				Rg2_z += (ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ())*(ingredients.getMolecules()[k].getZ()-ingredients.getMolecules()[l].getZ());
			}
			monomerCounter++;
		}
		if(monomerCounter != ingredients.getMolecules().size())
		{
			throw std::runtime_error("invalid number of monomers in Rg2-calculation");
		}

		Rg2_x /= double(ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_y /= double(ingredients.getMolecules().size()*ingredients.getMolecules().size());
		Rg2_z /= double(ingredients.getMolecules().size()*ingredients.getMolecules().size());

		Rg2 = Rg2_x+Rg2_y+Rg2_z;
		
        Statistic_Rg2.AddValue(Rg2); // radius of gyration
        Statistic_Rg2_x.AddValue(Rg2_x); // radius of gyration in x
        Statistic_Rg2_y.AddValue(Rg2_y); // radius of gyration in y
        Statistic_Rg2_z.AddValue(Rg2_z); // radius of gyration in z
        
    }
    
    return true;
}

    
struct PathSeparator {

    bool operator()(char ch) const {
        return ch == '\\' || ch == '/';
    }
};

template<class IngredientsType>
bool AnalyzerConstantForce<IngredientsType>::isFileExisting(std::string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

template<class IngredientsType>
void AnalyzerConstantForce<IngredientsType>::cleanup() {

    // print results into a file
    // get the filename and path
    // find the filename without path and extensions
    std::string filenameGeneral = std::string(std::find_if(ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator()).base(), ingredients.getName().end());

    std::string::size_type const p(filenameGeneral.find_last_of('.'));
    filenameGeneral = filenameGeneral.substr(0, p);

    std::stringstream systeminfo;
	systeminfo << "f0_" << std::setprecision (10) << std::fixed << ingredients.getAmplitudeForce() << "_";

    // construct a list
    std::vector < std::vector<double> > tmpResults_ConstantForce;

    // we have 41 columns and row
    uint32_t columns = 41;
    uint32_t rows = 1;



    // we have columns
    tmpResults_ConstantForce.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResults_ConstantForce[i].resize(rows);

   
    tmpResults_ConstantForce[0][0] = ingredients.getAmplitudeForce();
     
    tmpResults_ConstantForce[1][0] = Statistic_Rx.ReturnM1();
    tmpResults_ConstantForce[2][0] = Statistic_Rx.ReturnM2();
    tmpResults_ConstantForce[3][0] = Statistic_Rx.ReturnVar();
    tmpResults_ConstantForce[4][0] = 2.0 * Statistic_Rx.ReturnSigma() / std::sqrt(1.0 * Statistic_Rx.ReturnN());
    tmpResults_ConstantForce[5][0] = Statistic_Rx.ReturnN();

    tmpResults_ConstantForce[6][0] = Statistic_Ry.ReturnM1();
    tmpResults_ConstantForce[7][0] = Statistic_Ry.ReturnM2();
    tmpResults_ConstantForce[8][0] = Statistic_Ry.ReturnVar();
    tmpResults_ConstantForce[9][0] = 2.0 * Statistic_Ry.ReturnSigma() / std::sqrt(1.0 * Statistic_Ry.ReturnN());
    tmpResults_ConstantForce[10][0] = Statistic_Ry.ReturnN();

    tmpResults_ConstantForce[11][0] = Statistic_Rz.ReturnM1();
    tmpResults_ConstantForce[12][0] = Statistic_Rz.ReturnM2();
    tmpResults_ConstantForce[13][0] = Statistic_Rz.ReturnVar();
    tmpResults_ConstantForce[14][0] = 2.0 * Statistic_Rz.ReturnSigma() / std::sqrt(1.0 * Statistic_Rz.ReturnN());
    tmpResults_ConstantForce[15][0] = Statistic_Rz.ReturnN();
    
    tmpResults_ConstantForce[16][0] = Statistic_Ree.ReturnM1();
    tmpResults_ConstantForce[17][0] = Statistic_Ree.ReturnM2();
    tmpResults_ConstantForce[18][0] = Statistic_Ree.ReturnVar();
    tmpResults_ConstantForce[19][0] = 2.0 * Statistic_Ree.ReturnSigma() / std::sqrt(1.0 * Statistic_Ree.ReturnN());
    tmpResults_ConstantForce[20][0] = Statistic_Ree.ReturnN();

    tmpResults_ConstantForce[21][0] = Statistic_Rg2_x.ReturnM1();
    tmpResults_ConstantForce[22][0] = Statistic_Rg2_x.ReturnM2();
    tmpResults_ConstantForce[23][0] = Statistic_Rg2_x.ReturnVar();
    tmpResults_ConstantForce[24][0] = 2.0 * Statistic_Rg2_x.ReturnSigma() / std::sqrt(1.0 * Statistic_Rg2_x.ReturnN());
    tmpResults_ConstantForce[25][0] = Statistic_Rg2_x.ReturnN();

    tmpResults_ConstantForce[26][0] = Statistic_Rg2_y.ReturnM1();
    tmpResults_ConstantForce[27][0] = Statistic_Rg2_y.ReturnM2();
    tmpResults_ConstantForce[28][0] = Statistic_Rg2_y.ReturnVar();
    tmpResults_ConstantForce[29][0] = 2.0 * Statistic_Rg2_y.ReturnSigma() / std::sqrt(1.0 * Statistic_Rg2_y.ReturnN());
    tmpResults_ConstantForce[30][0] = Statistic_Rg2_y.ReturnN();

    tmpResults_ConstantForce[31][0] = Statistic_Rg2_z.ReturnM1();
    tmpResults_ConstantForce[32][0] = Statistic_Rg2_z.ReturnM2();
    tmpResults_ConstantForce[33][0] = Statistic_Rg2_z.ReturnVar();
    tmpResults_ConstantForce[34][0] = 2.0 * Statistic_Rg2_z.ReturnSigma() / std::sqrt(1.0 * Statistic_Rg2_z.ReturnN());
    tmpResults_ConstantForce[35][0] = Statistic_Rg2_z.ReturnN();

    tmpResults_ConstantForce[36][0] = Statistic_Rg2.ReturnM1();
    tmpResults_ConstantForce[37][0] = Statistic_Rg2.ReturnM2();
    tmpResults_ConstantForce[38][0] = Statistic_Rg2.ReturnVar();
    tmpResults_ConstantForce[39][0] = 2.0 * Statistic_Rg2.ReturnSigma() / std::sqrt(1.0 * Statistic_Rg2.ReturnN());
    tmpResults_ConstantForce[40][0] = Statistic_Rg2.ReturnN();



    

    std::stringstream comment;
    comment << " File produced by analyzer AnalyzerConstantForce" << std::endl
            << " for constant force in x-direction: f = f0 " << std::endl
            << std::endl
            << std::endl
            << " f0      ... offset force (col: 1)" << std::endl
            << " <|R|>   ... end-to-end-distance = <(R*R)^(0.5)>= <(Rx^2+Ry^2+Rz^2)^(0.5)> (col: 17-21)" << std::endl
            << " <Rx>    ... end-to-end-vector component in x-direction (col:  2-6)" << std::endl
            << " <Ry>    ... end-to-end-vector component in y-direction (col:  7-11)" << std::endl
            << " <Rz>    ... end-to-end-vector component in z-direction (col: 12-16)" << std::endl
            << " <Rg2>   ... radius of gyration = 1/N sum (vec r_i - vec r_COM)^2 (col: 37-41)" << std::endl
            << " <Rg2_x> ... radius of gyration component in x-direction (col: 22-26)" << std::endl
            << " <Rg2_y> ... radius of gyration component in y-direction (col: 27-31)" << std::endl
            << " <Rg2_z> ... radius of gyration component in z-direction (col: 32-36)" << std::endl
            << " under N = " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << std::endl
            << "f0 <Rx> <Rx²> Varianz(Rx) Error(Rx) SampleSize(Rx) <Ry> <Ry²> Varianz(Ry) Error(Ry) SampleSize(Ry) <Rz> <Rz²> Varianz(Rz) Error(Rz) SampleSize(Rz) <|R|> <|R|²> Varianz(|R|) Error(|R|) SampleSize(|R|) <Rg2_x> <(Rg2_x)²> Varianz(Rg2_x) Error(Rg2_x) SampleSize(Rg2_x) <Rg2_y> <(Rg2_y)²> Varianz(Rg2_y) Error(Rg2_y) SampleSize(Rg2_y) <Rg2_z> <(Rg2_z)²> Varianz(Rg2_z) Error(Rg2_z) SampleSize(Rg2_z) <Rg2> <(Rg2)²> Varianz(Rg2) Error(Rg2) SampleSize(Rg2)";
            


    //new filename
    std::string filename_ConstantForceReeRg2 = filenameGeneral + "_" + systeminfo.str() + "_ReeRg2.dat";
    
    if(!isFileExisting(dstdir+"/"+filename_ConstantForceReeRg2))
        ResultFormattingTools::writeResultFile(dstdir+"/"+filename_ConstantForceReeRg2, this->ingredients, tmpResults_ConstantForce, comment.str());
    else
        ResultFormattingTools::appendToResultFile(dstdir + "/" + filename_ConstantForceReeRg2, tmpResults_ConstantForce);

    
    // ====================================================
    // List all entries of hysteresis area
    // ====================================================
    
    std::vector < std::vector<double> > tmpResults_ConstantForceIntegralList;

    // we have columns and row
    columns = 1;
    rows = vectorAllRxx.size();

    // we have columns
    tmpResults_ConstantForceIntegralList.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResults_ConstantForceIntegralList[i].resize(rows);
    
    // C++11 alternative:
    size_t counter_vecReeList = 0;
    
    for (const auto& entry : vectorAllRxx) {
       
        tmpResults_ConstantForceIntegralList[0][counter_vecReeList] = entry;
        
        counter_vecReeList++;
    }

    std::stringstream comment_ReeList;
    comment_ReeList << " File produced by analyzer AnalyzerConstantForce" << std::endl
            << " for constant force in x-direction: f = f0 " << std::endl
            << std::endl
            << std::endl
            << " f0    ... offset force (col: 1)" << std::endl
            << " Rx    ... end-to-end-vector component in x-direction (col: 1)" << std::endl
            << " under N = " << ingredients.getMolecules().size() << " monomers" << std::endl
            << std::endl
            << std::endl
            << "Rx";
           

    //new filename
    std::string filename_ConstantForceReeRg2Integral = filenameGeneral + "_" + systeminfo.str() + "_ReeList.dat";
    
    if(!isFileExisting(dstdir+"/"+filename_ConstantForceReeRg2Integral))
        ResultFormattingTools::writeResultFile(dstdir+"/"+filename_ConstantForceReeRg2Integral, this->ingredients, tmpResults_ConstantForceIntegralList, comment_ReeList.str());
    else
        ResultFormattingTools::appendToResultFile(dstdir + "/" + filename_ConstantForceReeRg2Integral, tmpResults_ConstantForceIntegralList);
    
    

}

#endif /*ANALYZER_CONSTANT_FORCE_H*/
