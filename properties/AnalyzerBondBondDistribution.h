
// HEADER
#ifndef ANALYZER_BOND_BOND_DISTRIBUTION_H
#define ANALYZER_BOND_BOND_DISTRIBUTION_H

// INCLUDE
#include <string>
#include <iostream>
#include <cmath>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DistanceCalculation.h>

// not in LeMonADE package
#include <HistogramGeneralStatistik1D.h>

template < class IngredientsType > class AnalyzerBondBondDistribution : public AbstractAnalyzer
{

private:
    //! lower histrogram border
    double low_dist_bound = 0.;
    //! upper histogram border
    double upp_dist_bound = 2. * M_PI;
    //! number of histrgarm bins
    int num_bins = 50;
    //! typedef for the underlying container holding the monomers
    typedef typename IngredientsType::molecules_type molecules_type;
    //! reference to the complete system
    const IngredientsType& ingredients;
    //! bb is calculated for the groups in this vector
    std::vector<MonomerGroup<molecules_type> > groups;
    //! histogram that accummulates the bond bond angle values
    HistogramGeneralStatistik1D bbdist;
    //! name of output files are outputFilePrefix_averages.dat and outputFilePrefix_timeseries.dat
    std::string outputFile;
    //! calculate the bond bond distrubution  of the monomer group
    void cummulateBBD(const MonomerGroup<molecules_type>& group) const;
    //! analyze only after equilibration time has been reached
    uint32_t equilibrationTime;
protected:
    //! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
    void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

    //! constructor
    AnalyzerBondBondDistribution(const IngredientsType& ing,
                             std::string filename="BondBondDistribution.dat",
                             uint32_t equilibrationTime_=0);

    //only used to make sure you initialize your groups before you do things
    bool initialized;
    //! destructor. does nothing
    virtual ~AnalyzerBondBondDistribution(){}
    //! Initializes data structures, empty histogram n bins. Called by TaskManager::initialize()
    virtual void initialize();
    //! calculates bond angles between all monomers. Called by TaskManager::execute()
    virtual bool execute();
    //! Writes the final results to file
    virtual void cleanup();
    //! Change the output file name
    void setOutputFile(std::string filename){outputFile=filename;isFirstFileDump=true;}
    //! Set equilibration time
    void setEquilibrationTime(uint32_t time){equilibrationTime=time;}
    //! Get equilibration time
    uint32_t getEquilibrationTime(){return equilibrationTime;}

};

// constructor
template<class IngredientsType>
AnalyzerBondBondDistribution<IngredientsType>::AnalyzerBondBondDistribution(const IngredientsType& ing, std::string outputFile_, uint64_t equilibrationTime_)
: ingredients(ing), outputFile(outputFIle_), equilibrationTime(equilibrationTime_)
{initialized=false;}

// initlaizer
template<class IngredientsType>
void AnalyzerBondBondDistribution<IngredientsType>::initialize() {

    // initialize histogram

    // ready to analyze
    initialized=true;
}

// exectue
// loop through monomer groups and calculate the bondbond distriubition
template<class IngredientsType>
bool AnalyzerBondBondDistribution<IngredientsType>::execute() {

    //check if groups have been initialized. if not, exit and explain
    if(initialized==false)
    {
        std::stringstream errormessage;
        errormessage<<"AnalyzerBondBondDistribution::execute()...analyzer not initialized\n"
        <<"Use AnalyzerBondBondDistribution::initialize() or Taskmanager::init()\n";

        throw std::runtime_error(errormessage.str());
    }

    // if the equilibriation time has been reached, evalulate
    if(ingredients.getMolecules().getAge() >= equilibrationTime)
    {
        // loop through each monomer group (i.e. polymer)
        // calculate the bond-bond angle for each neighboring bond-bond pair
        // accumulate in the histogram

    }

    return true;
}

template<class IngredientsType>
void AnalyzerBondBondDistribution<IngredientsType>::cleanup() {

    // print results into a file

}


