
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
    //! typedef for the underlying container holding the monomers
    typedef typename IngredientsType::molecules_type molecules_type;
    //! reference to the complete system
    const IngredientsType& ingredients;
    //! bb is calculated for the groups in this vector
    std::vector<MonomerGroup<molecules_type> > groups;
    //! histogram that accummulates the bond bond angle values
    HistogramGeneralStatistik1D
    //! timeseries of the Re2e. Components: [0]-> Re2e_x, [1]->Re2e_y, [2]->Re2e_z, [3]:Re2e_tot
    std::vector< std::vector<double> > Re2eTimeSeries;
    //! vector of mcs times for writing the time series
    std::vector<double> MCSTimes;
    //! max length of internal buffer for each coordinate, before saving to disk
    uint32_t bufferSize;
    //! name of output files are outputFilePrefix_averages.dat and outputFilePrefix_timeseries.dat
    std::string outputFile;
    //! flag used in dumping time series output
    bool isFirstFileDump;
    //! save the current values in Rg2TimeSeriesX, etc., to disk
    void dumpTimeSeries();
    //! calculate the Rg squared of the monomer group
    VectorDouble3 calculateRe2eComponents(const MonomerGroup<molecules_type>& group) const;
    //! starting number of MCS when analyzer is initialized
    uint32_t startMCS;
    //! analyze after equilibration time
    uint32_t equilibrationTime;
protected:
    //! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
    void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

    //! constructor
    AnalyzerBondBondDistribution(const IngredientsType& ing,
                             std::string filename="BondBondDistribution.dat",
                             uint32_t equilibrationTime_=0);

    //! destructor. does nothing
    virtual ~AnalyzerBondBondDistribution(){}
    //! Initializes data structures, empty histogram n bins. Called by TaskManager::initialize()
    virtual void initialize();
    //! calculates bond angles between all monomers. Called by TaskManager::execute()
    virtual bool execute();
    //! Writes the final results to file
    virtual void cleanup();
    //! Set the number of values, after which the time series is saved to disk
    void setBufferSize(uint32_t size){bufferSize=size;}
    //! Change the output file name
    void setOutputFile(std::string filename){outputFile=filename;isFirstFileDump=true;}
    //! Set equilibration time
    void setEquilibrationTime(uint32_t time){equilibrationTime=time;}
    //! Get equilibration time
    uint32_t getEquilibrationTime(){return equilibrationTime;}

};


