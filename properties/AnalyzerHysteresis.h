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

// from LeoMonADE
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>

// not in LeMonADE package
#include <HistogramGeneralStatistik1D.h>

/* **************************************************************
 * A simple analyzer example for static properties 
 * 
 * In this case it calculates the average RgSquared of the molecules in
 * the system. Plain and simple...well...at least plain. Have fun.
 * *************************************************************/

template<class IngredientsType>
class AnalyzerHysteresis : public AbstractAnalyzer {
public:
    // default constructor
    AnalyzerHysteresis(const IngredientsType& ing, std::string filename_, uint64_t evalulation_time_, uint64_t save_interval_);
    // deconstructor
    virtual ~AnalyzerHysteresis() {};
    //typedef typename IngredientsType::molecules_type molecules_type;
    const typename IngredientsType::molecules_type& molecules;
    //! initialize analyzer
    virtual void initialize();
    //! execute analyzer
    virtual bool execute();
    //! finish analysis, write to file
    virtual void cleanup();
    // check if file exists
    bool isFileExisting(std::string);

private:
    // constant tolerance used for bounds
    const double TOL = 0.001;
    // constant, maximum possible value from hysteresis end-to-end distance
    const double hys_hist_max_value = 1.;
    // constant, minimum possible value from hysteresis end-to-end distance
    const double hys_hist_min_value = -1.;
    // constant, number of bins in hysteresis histogram (between min and max values)
    const int hys_hist_bins = 200.;
    // length used to normalize end-to-end distance, which is the total chain length
    double norm_length;
    // counts the number of times that execution has iterated
    int hys_counter;
    // counts the number of times that a hysteresis loop is completed
    int loop_counter;
    // points to LeMonADE ingredients
    const IngredientsType& ingredients;
    //only used to make sure you initialize your groups before you do things
    bool initialized;
    // used to determine if analyzer is ready to start collecting data
    bool equilibriated;
    //! Time in MCS, which need to pass to calculate the observables
    uint64_t equilibriation_time;
    //! number of points that are sampled along the hysteresis loop
    int numPeriodPoints;
    //! Time in MCS within one period, there the observables is collected
    uint64_t step_interval;
    //! period of hysteresis loop, in MCSs
    int hys_period;
    //! vector that total distance is projected onto, cooresponds to orientation of force
    VectorDouble3 projVec;
    //! amplitude of oscillatory force
    double hys_amplitude;
    //! intermediate array that collects the average end to end distance over one period
    std::vector<double> single_hys_loop;
    //! array of histogram statistics for each point sampled along hysteresis curve over course of simulation
    std::vector<HistogramGeneralStatistik1D> averaged_hys_loop;
    //! collects statistics for hysteresis integral, averaged for each loop
    std::vector<double> avg_hys_integral;
    // monomers experience forcing
    uint32_t indexMonomersForce[2];
    // calculation of hysteresis within on period using the trapezoidal method
    double calculateHysteresisInOnePeriodWithTrapezoidalMethod();
    // outfile name
    std::string outfile;
    // maximum hysteresis value calculated over the course of the simulation loop
    double max_hys_val;
    // minimum hysteresis value calculated over the course of the simulation loop
    double min_hys_val;

    uint32_t counterFrames;
    std::string filename;
    std::string dstdir;

    
};




/*****************************************************************************
 * IMPLEMENTATION OF METHODS
 * **************************************************************************/


/////////////////////////////////////////////////////////////////////////////

template<class IngredientsType>
AnalyzerHysteresis<IngredientsType>::AnalyzerHysteresis(const IngredientsType& ing, std::string filename_, uint64_t evalulation_time_, uint64_t save_interval_)
: ingredients(ing), molecules(ing.getMolecules()), outfile(filename_), equilibriation_time(evalulation_time_), step_interval(save_interval_) {
    counterFrames = 1;
    initialized=false;
    equilibriated = false;
    hys_counter = 0;
    loop_counter = 0;
    min_hys_val = 0;
    max_hys_val = 0;
}

template<class IngredientsType>
void AnalyzerHysteresis<IngredientsType>::initialize() {

    // check that the force oscillation is on
    if (! ingredients.isForceOscillationOn()){
        // if the force oscillation is not on, throw an error
        // the hysteresis analyzer cannot run if the force oscillation is not on
        std::stringstream errormessage;
        errormessage << "AnalyzerHysteresis::initialize()...Force Oscillation is off. Hysteresis analyzer cannot work." << std::endl;
        throw std::runtime_error(errormessage.str());
    }

    // check that the period is alligned with the sampling frequency
    hys_period = ingredients.getForceOscillationPeriod();
    hys_amplitude = ingredients.getForceOscillationAmplitude();
    projVec = ingredients.getForceVector();
    if (hys_period % step_interval != 0) {
        // if the period of the save interval does not evenly go into the hysteresis period
        // throw an error because that will cause some issues in the fututre
        std::stringstream errormessage;
        errormessage << "AnalyzerHysteresis::initialize()...Force Oscillation Period (" << hys_period << " MCS) and save interval period (" << step_interval << " MCS) not evenly divisable. Integer representing number of save points along hysteresis loop needs to be an even number." << std::endl;
        throw std::runtime_error(errormessage.str());
    } else {
        // the numbers are evenly divisable
        numPeriodPoints = hys_period / step_interval;
    }

    //get the monomers with attributes subjected to force
	uint32_t n = 0;
	for (uint32_t i=0; i < ingredients.getMolecules().size();i++){
		if (ingredients.getMolecules()[i].getAttributeTag() != 0){
			indexMonomersForce[n] = i;
			n += 1;
		}
    }
    // only two monomers can be subjected to a force
    if((n > 2) || (n==0))
    {
        std::stringstream errormessage;
        errormessage<<"AnalyzerHysteresis::initialize()...ambigious number of monomers subjected to force" << std::endl;

        throw std::runtime_error(errormessage.str());
    }

    // determine if the polymer is a ring or a chain .. (or a chain of rings .. ?)
    int num_monomers = ingredients.getMolecules().size();
    if (indexMonomersForce[0] == 0) {
        // if the first momonmer experiencing force is the first in the chain, check the second
        // NOTE :: this should always be the case
        if (indexMonomersForce[1] == num_monomers - 1) {
            // the polymer is a linear chain, the maximum length is the total number of momonmers
            norm_length = (double) num_monomers;
        } else {
            // the polymer is a ring, the maximum length is half the total number of monomers
            norm_length = (double) num_monomers / 2.;
        }
    }
    norm_length = 3. * norm_length;

    // initialize the histograms / arrays that collect hysteresis statistics
    single_hys_loop.resize(numPeriodPoints);
    averaged_hys_loop.resize(numPeriodPoints);
    for (int n; n < averaged_hys_loop.size(); n++) {
        // initialize each of the arrays for each of the points along the period that will be sampled
        // one histogram array represents the averaged end-to-end distance at each time point along the hystersis (oscillation) loop
        averaged_hys_loop[n] = HistogramGeneralStatistik1D(hys_hist_min_value, hys_hist_max_value, hys_hist_bins);
        // one more histogram represents hysteresis which is calculated after the completion of each loop
        single_hys_loop[n] = 0.;
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

    if (not equilibriated) {
        // check if the system is ready to start tracking
        int tnow = ingredients.getMolecules().getAge();
        equilibriated = (tnow >= equilibriation_time) && ((tnow % hys_period) == 0);
    } else {

        // iterate hysteresis counter
        hys_counter++;
        int period_idx = (hys_counter - 1) % numPeriodPoints; // index corresponding to the current phase in period
        
        // calculate the end-to-end vector of the polymer chain / ring projected onto the force vector
        VectorDouble3 diff_vec;
        diff_vec.setX((ingredients.getMolecules()[indexMonomersForce[1]]).getX() - (ingredients.getMolecules()[indexMonomersForce[0]]).getX());
        diff_vec.setY((ingredients.getMolecules()[indexMonomersForce[1]]).getY() - (ingredients.getMolecules()[indexMonomersForce[0]]).getY());
        diff_vec.setZ((ingredients.getMolecules()[indexMonomersForce[1]]).getZ() - (ingredients.getMolecules()[indexMonomersForce[0]]).getZ());
        double mag = (diff_vec * projVec) / (projVec * projVec);
        diff_vec = mag * projVec;

        // calculate the normalized distance
        double len = std::sqrt(std::pow(diff_vec.getX(), 2)+std::pow(diff_vec.getY(),2)+std::pow(diff_vec.getZ(),2)) / norm_length;
        mag = diff_vec * projVec;
        if (mag < 0) {
            // if the difference vector is negative relative to the force vector, the length is negative
            len = - len;
        }

        // accumulate the length according to the current phase of oscillation
        single_hys_loop[period_idx] = len;

        // if a complete period cycle has been reached
        if (period_idx == (numPeriodPoints - 1)) {
            // iterate the hysteresis loop counter
            loop_counter++;
            // calculate hystresis, accumulate
            double hysteresis_integral = calculateHysteresisInOnePeriodWithTrapezoidalMethod();
            avg_hys_integral.push_back(hysteresis_integral);
            // avg_single_hys_loop.addValue(hysteresis_integral);

            // check min and max values
            if (loop_counter == 1) {
                // if first loop, initialize min and max values
                min_hys_val = hysteresis_integral;
                max_hys_val = hysteresis_integral;
            } else {
                // check for a new maximum or minimum
                if (hysteresis_integral > max_hys_val) {max_hys_val = hysteresis_integral;}
                if (hysteresis_integral < min_hys_val) {min_hys_val = hysteresis_integral;}
            }

            // accumulate the loop in the histogram for hysteresis, reset and repeat
            for (int n; n < averaged_hys_loop.size(); n++) {
                averaged_hys_loop[n].addValue(single_hys_loop[n], 1.); // add the value to the histogram tracking statistics around the hysteresis loop
                single_hys_loop[n] = 0.; // reset for the next loop
            }
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

    double omega = 2.0 * 3.14159265358979323846 / hys_period;
    // double fA = ingredients.getAmplitudeOscillatoryForce();
    /* Finding Integration Value */

    // f = f0+fA*sin(omega*t) in eX-direction

    //         /         /T
    //        |          |
    // A(T) = o R * df = o R(f)*fA*omega*cos(omega*t) dt
    //        |          |
    //        /         0/
 
    double hysteresis_integral = 0.0;
 
    for(size_t i = 0; i < single_hys_loop.size(); i++) {
        int i_1 = i;
        int i_2 = i + 1;
        if (i_2 == numPeriodPoints) {i_2 = 0;}
        // double t_1 = i_1 * stepSize;
        // double t_2 = i_2 * stepSize;
        double F_1 = ingredients.getForceNow(i_1 * stepSize);
        double F_2 = ingredients.getForceNow(i_2 * stepSize);
        // hysteresis_integral += single_hys_loop.at(i) * omega * hys_amplitude * std::cos(omega * i * stepSize) * stepSize;
        hysteresis_integral += 0.5 * (single_hys_loop[i_2] - single_hys_loop[i_1]) * (F_2 - F_1);
    }
 
    // The first value counts twice (2*0.5=1)at this is Rxx at t%T=0 in the cycle
    // cos(omega * 0 * stepSize) = cos(omega * T * stepSize) != 1
    // hysteresis_integral += AllRxxWithinPeriod.at(0) * omega * hys_amplitude * stepSize;

    return hysteresis_integral;
}
    
struct PathSeparator
{
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

    //  create histogram contining all hystresis loop values
    HistogramGeneralStatistik1D avg_single_hys_loop = HistogramGeneralStatistik1D(min_hys_val, max_hys_val, hys_hist_bins);
    for (int n = 0; n < avg_hys_integral.size(); n++) {
        avg_single_hys_loop.addValue(avg_hys_integral[n], 1.);
    }

    // if the file exists, delete it
    char* outfile_char_array = new char[outfile.length() + 1];
    strcpy(outfile_char_array, outfile.c_str());
    remove(outfile_char_array);

    // open file stream
    std::ofstream file;
    file.open(outfile.c_str());

    // write the header
    file << "i,f,avg,std,n" << std::endl;

    // print the average hysteresis value
    file << (0) << ",NA," << avg_single_hys_loop.getHistAverage() << ",NA," << avg_single_hys_loop.getNCounts() << std::endl;

    // print results into a file
    for (int n = 0; n < averaged_hys_loop.size(); n++) {
        file << (n+1) << "," << ingredients.getForceNow(n * step_interval) << "," << averaged_hys_loop[n].getHistAverage() << ",NA," << averaged_hys_loop[n].getNCounts() << std::endl;
    }
    file.close();

    // print results into a file
    // get the filename and path
    // find the filename without path and extensions
    std::string filenameGeneral = std::string(std::find_if(ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator()).base(), ingredients.getName().end());
/*
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
        ResultFormattingTools::appendToResultFile(dstdir + "/" + filename_HysteresisIntegral, tmpResults_HysteresisIntegralList);*/
    
    

}

#endif /*ANALYZERHYSTERESIS_H*/
