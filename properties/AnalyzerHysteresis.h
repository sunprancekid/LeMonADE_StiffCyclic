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
    uint64_t hys_period;
    //! vector that total distance is projected onto, cooresponds to orientation of force
    VectorDouble3 projVec;
    //! amplitude of oscillatory force
    double hys_amplitude;
    //! intermediate array that collects the measured end to end distance over one period
    std::vector<double> single_hys_loop_distance;
    //! intermediate array that collects the current force that a chain is experiencing
    std::vector<double> single_hys_loop_force;
    //! array of containing the averaged end to end distance for each point sampled along period over course of simulation
    std::vector<HistogramGeneralStatistik1D> avg_hys_loop_dist;
    //! array of containing the averaged force for each point sampled along period over course of simulation
    std::vector<HistogramGeneralStatistik1D> avg_hys_loop_forc;
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
    // norm_length = 1.;

    // initialize the histograms / arrays that collect hysteresis statistics
    // arrays that collect distance / force over one loop
    single_hys_loop_distance.resize(numPeriodPoints);
    single_hys_loop_force.resize(numPeriodPoints);
    // array that contains averaged end to end distance along the hystersis curve
    avg_hys_loop_dist.resize(numPeriodPoints);
    for (int n; n < avg_hys_loop_dist.size(); n++) {
        // initialize each of the arrays for each of the points along the period that will be sampled
        // one histogram array represents the averaged end-to-end distance at each time point along the hystersis (oscillation) loop
        avg_hys_loop_dist[n] = HistogramGeneralStatistik1D(hys_hist_min_value, hys_hist_max_value, hys_hist_bins);
        // one more histogram represents hysteresis which is calculated after the completion of each loop
        // single_hys_loop[n] = 0.;
    }
    // array that contains averaged force experience by chain along the hysteresis curve
    avg_hys_loop_forc.resize(numPeriodPoints);
    for (int n; n < avg_hys_loop_forc.size(); n++) {
        // initilize each histogram with the minimum and maximum forces that the chain can experience
        // (determined by amplitude)
        avg_hys_loop_forc[n] = HistogramGeneralStatistik1D(-hys_amplitude - TOL, hys_amplitude + TOL, hys_hist_bins);
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
        // in order to start tracking, the sysmtem must be at the beginning of the period
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
        double len = std::sqrt(std::pow(diff_vec.getX(), 2) + std::pow(diff_vec.getY(), 2) + std::pow(diff_vec.getZ(), 2)) / norm_length;
        mag = diff_vec * projVec;
        if (mag < 0) {
            // if the difference vector is negative relative to the force vector, the length is negative
            len = - len;
        }

        // accumulate the length according to the current phase of oscillation
        single_hys_loop_distance[period_idx] = len; // actual end-to-end distance relative to force vector
        single_hys_loop_force[period_idx] = ingredients.getForceNow(ingredients.getMolecules().getAge() % hys_period);
        // get the actual force experienced by polymer, use mod to prevent errors with over counting

        // if a complete period cycle has been reached
        if (period_idx == (numPeriodPoints - 1)) {

            // iterate the hysteresis loop counter
            loop_counter++;

            // calculate hystresis, accumulate
            double hysteresis_integral = calculateHysteresisInOnePeriodWithTrapezoidalMethod();
            avg_hys_integral.push_back(hysteresis_integral); // add calculated hysteresis
            // avg_single_hys_loop.addValue(hysteresis_integral);

            // TODO :: if analyzer works properly, hystersis value should always be between zero and one
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

            // accumulate the average end to end distance / force along the hysteresis loop, reset and repeat
            for (int n; n < avg_hys_loop_dist.size(); n++) {
                // add the value to the histogram tracking statistics around the hysteresis loop
                avg_hys_loop_dist[n].addValue(single_hys_loop_distance[n], 1.);
                single_hys_loop_distance[n] = 0.; // reset for the next loop
                // collect the force along the loop length
                avg_hys_loop_forc[n].addValue(single_hys_loop_force[n], 1.);
                single_hys_loop_force[n] = 0.; // reset for next loop
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
 
    for(size_t i = 0; i < single_hys_loop_distance.size(); i++) {
        // establish integers that access arrays
        int i_1 = i;
        int i_2 = i + 1;
        if (i_2 == numPeriodPoints) {i_2 = 0;}
        // calculate hysteresis according to trapazoidal method
        hysteresis_integral += 0.5 * (single_hys_loop_distance[i_2] + single_hys_loop_distance[i_1]) * (single_hys_loop_force[i_2] - single_hys_loop_force[i_1]);
    }

    hsyteresis_integral = hys_amplitude * omega * hsyteresis_integral;
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
    file << "i,f_exp,f_avg,f_std,E2E_avg,E2E_std,n" << std::endl;

    // print the average hysteresis value
    file << (0) << ",NA,NA,NA," << avg_single_hys_loop.getHistAverage() << "," << avg_single_hys_loop.getHistSTD() << "," << avg_single_hys_loop.getNCounts() << std::endl;

    // print results into a file
    for (int n = 0; n < avg_hys_loop_dist.size(); n++) {
        file << (n+1) << "," << ingredients.getForceNow(n * step_interval) << "," << avg_hys_loop_forc[n].getHistAverage() << "," << avg_hys_loop_forc[n].getHistSTD() << "," << avg_hys_loop_dist[n].getHistAverage() << "," << avg_hys_loop_dist[n].getHistSTD() << "," << avg_hys_loop_dist[n].getNCounts() << std::endl;
    }
    file.close();

    // print results into a file
    // get the filename and path
    // find the filename without path and extensions
    std::string filenameGeneral = std::string(std::find_if(ingredients.getName().rbegin(), ingredients.getName().rend(), PathSeparator()).base(), ingredients.getName().end());
}

#endif /*ANALYZERHYSTERESIS_H*/
