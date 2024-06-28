
// HEADER
#ifndef ANALYZER_BOND_BOND_DISTRIBUTION_H
#define ANALYZER_BOND_BOND_DISTRIBUTION_H

// INCLUDE
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
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
    const double low_dist_bound = -1.;
    //! upper histogram border
    const double upp_dist_bound = 1.;
    //! actual number of histogram bins
    int num_bins;
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
    std::vector<double> cummulateBBD(const MonomerGroup<molecules_type>& group) const;
    //! analyze only after equilibration time has been reached
    uint32_t equilibrationTime;
protected:
    //! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
    void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

    //! constructor
    AnalyzerBondBondDistribution(const IngredientsType& ing,
                                 std::string filename="BondBondDistribution.dat",
                                 uint32_t equilibrationTime_=0,
                                 int nbins_ = 25);
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
    void setOutputFile(std::string filename){outputFile=filename;}
    //! Set equilibration time
    void setEquilibrationTime(uint32_t time){equilibrationTime=time;}
    //! Get equilibration time
    uint32_t getEquilibrationTime(){return equilibrationTime;}

};

// constructor
template<class IngredientsType>
AnalyzerBondBondDistribution<IngredientsType>::AnalyzerBondBondDistribution(const IngredientsType& ing, std::string filename, uint32_t equilibrationTime_, int nbins_)
: ingredients(ing), outputFile(filename), equilibrationTime(equilibrationTime_), num_bins(nbins_)
{initialized=false;}

// initlaizer
template<class IngredientsType>
void AnalyzerBondBondDistribution<IngredientsType>::initialize()
{

    // initialize histogram
    bbdist = HistogramGeneralStatistik1D(low_dist_bound, upp_dist_bound, num_bins);
    //if no groups are set, use the complete system by default
    //groups can be set using the provided access function
    if(groups.size()==0)
    {
        groups.push_back(MonomerGroup<molecules_type>(ingredients.getMolecules()));
        for(size_t n=0;n<ingredients.getMolecules().size();n++)
            groups[0].push_back(n);
    }
    // ready to analyze
    initialized=true;
}

// exectue
// loop through monomer groups and calculate the bondbond distriubition
template<class IngredientsType>
bool AnalyzerBondBondDistribution<IngredientsType>::execute()
{

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
        // loop through each monomer group (i.e. polymer)for(size_t n=0;n<groups.size();n++)
        for(size_t n=0;n<groups.size();n++) {
            // calculate the bond-bond angle for each neighboring bond-bond pair
            // accumulate in the histogram
            std::vector<double> angles = cummulateBBD(groups[n]);
            for (int m = 0; m < angles.size(); m++) {
                bbdist.addValue(angles[m], 1.);
            }
        }
    }

    return true;
}

template<class IngredientsType>
void AnalyzerBondBondDistribution<IngredientsType>::cleanup()
{

    // open file stream
    std::ofstream file;
    file.open(outputFile.c_str());

    // TODO check if file is open, etc.

    // print results into a file
    double sum = 0.;
    for (int n = 0; n < bbdist.getNBins(); n++) {
        sum += bbdist.getFirstMomentInBin(n);
        file << bbdist.getCenterOfBin(n) << ", " << bbdist.getFirstMomentInBin(n) << std::endl;
    }
    file.close();
}

// cummulate bond bond angles for polymer chain
template<class IngredientsType>
std::vector<double> AnalyzerBondBondDistribution<IngredientsType>::cummulateBBD(
    const MonomerGroup<molecules_type>& group) const
    {
        // get the group size
        double group_size = group.size();

        // list of angles
        std::vector<double> angles;

        //if group is empty, return zero vector
        if(group_size==0){
            return angles;
        }

        angles.resize(group_size - 1);
        // loop through each neighboring bond-bond pair
        for (int n = 0; n < (group_size - 2); n++){
            // get the vectors representing the first and second bonds
            VectorDouble3 bi;
            bi.setX(group[n+1].getX() - group[n].getX());
            bi.setY(group[n+1].getY() - group[n].getY());
            bi.setZ(group[n+1].getZ() - group[n].getZ());
            VectorDouble3 bj;
            bj.setX(group[n+2].getX() - group[n+1].getX());
            bj.setY(group[n+2].getY() - group[n+1].getY());
            bj.setZ(group[n+2].getZ() - group[n+1].getZ());
            // calculate the angle between the two bonds
            double angle;
            angle = bi*bi;
            angle *= bj*bj;
            angle = sqrt(angle);
            angle = (bi*bj)/angle;
            angles[n] = angle;
        }
    return angles;
    }

#endif /*ANALYZER_BOND_BOND_DISTRIBUTION_H*/


