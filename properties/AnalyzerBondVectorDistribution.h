
// HEADER
#ifndef ANALYZER_BOND_VECTOR_DISTRIBUTION_H
#define ANALYZER_BOND_VECTOR_DISTRIBUTION_H

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
#include <LeMonADE/utility/FastBondset.h>

// not in LeMonADE package
#include <HistogramGeneralStatistik1D.h>

template < class IngredientsType > class AnalyzerBondVectorDistribution : public AbstractAnalyzer
{

private:
    //! lower histrogram border
    const double low_dist_bound = -0.5;
    //! upper histogram border
    const double upp_dist_bound = 5.5;
    //! actual number of histogram bins
    const int num_bins = 6;
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
    std::vector<double> cummulateBVD(const MonomerGroup<molecules_type>& group) const;
    //! analyze only after equilibration time has been reached
    uint32_t equilibrationTime;
    //! contains bondset vectors
    FastBondset fb;
protected:
    //! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
    void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

    //! constructor
    AnalyzerBondVectorDistribution(const IngredientsType& ing,
                                 std::string filename="BondVectorDistribution.dat",
                                 uint32_t equilibrationTime_=0);
    //only used to make sure you initialize your groups before you do things
    bool initialized;
    //! destructor. does nothing
    virtual ~AnalyzerBondVectorDistribution(){}
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
AnalyzerBondVectorDistribution<IngredientsType>::AnalyzerBondVectorDistribution(const IngredientsType& ing, std::string filename, uint32_t equilibrationTime_)
: ingredients(ing), outputFile(filename), equilibrationTime(equilibrationTime_)
{initialized=false;}

// initlaizer
template<class IngredientsType>
void AnalyzerBondVectorDistribution<IngredientsType>::initialize()
{
    // initialize bondset vectors
    fb = FastBondset();
    fb.addBFMclassicBondset();
    fb.updateLookupTable();
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
bool AnalyzerBondVectorDistribution<IngredientsType>::execute()
{

    //check if groups have been initialized. if not, exit and explain
    if(initialized==false)
    {
        std::stringstream errormessage;
        errormessage<<"AnalyzerBondVectorDistribution::execute()...analyzer not initialized\n"
        <<"Use AnalyzerBondVectorDistribution::initialize() or Taskmanager::init()\n";

        throw std::runtime_error(errormessage.str());
    }

    // if the equilibriation time has been reached, evalulate
    if(ingredients.getMolecules().getAge() >= equilibrationTime)
    {
        // loop through each monomer group (i.e. polymer)for(size_t n=0;n<groups.size();n++)
        for(size_t n=0;n<groups.size();n++) {
            // calculate the bond-bond angle for each neighboring bond-bond pair
            // accumulate in the histogram
            std::vector<double> angles = cummulateBVD(groups[n]);
            for (int m = 0; m < angles.size(); m++) {
                bbdist.addValue(angles[m], 1.);
            }
        }
    }

    return true;
}

template<class IngredientsType>
void AnalyzerBondVectorDistribution<IngredientsType>::cleanup()
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
std::vector<double> AnalyzerBondVectorDistribution<IngredientsType>::cummulateBVD(
    const MonomerGroup<molecules_type>& group) const
    {
        // get the group size
        double group_size = group.size();

        // list of bond vectors catagorized by symmetry
        std::vector<double> bondvecs;

        //if group is empty, return zero vector
        if(group_size==0) {
            return bondvecs;
        }

        bondvecs.resize(group_size - 1);
        // loop through each neighboring bond-bond pair
        for (int n = 0; n < (group_size - 2); n++){
            // get the vector representing the bond
            VectorDouble3 bv;
            bv.setX(group[n+1].getX() - group[n].getX());
            bv.setY(group[n+1].getY() - group[n].getY());
            bv.setZ(group[n+1].getZ() - group[n].getZ());
            // translate the bond vector to an ASCII key
            int key = fb.getBondIdentifier(bv.getX(), bv.getY(), bv.getZ());
            // accumulate the vector according to the key
            // TODO merge (3 0 0) and (2 0 0), b/c same orientation
            if (key < 23) {
                // P +- (2 0 0)
                bondvecs[n] = 0;
            } else if (key < 47) {
                // P +- (1 2 0)
                bondvecs[n] = 1;
            } else if (key < 71) {
                // P +- (2 1 1)
                bondvecs[n] = 2;
            } else if (key < 77) {
                // P +- (3 0 0)
                bondvecs[n] = 3;
            } else if (key < 101) {
                // P +- (2 2 1)
                bondvecs[n] = 4;
            } else if (key < 125) {
                // P +- (3 1 0)
                bondvecs[n] = 5;
            } else{
                std::cout << "ERROR :: BondVectorDistribution :: key (" << key << ") outside of bondvector range." << std::endl;
                exit(0);
            }
        }
    return bondvecs;
    }

#endif /*ANALYZER_BOND_VECTOR_DISTRIBUTION_H*/


