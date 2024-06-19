
// HEADER
#ifndef ANALYZER_BOND_BOND_CORRELATION_H
#define ANALYZER_BOND_BOND_CORRELATION_H

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

template < class IngredientsType > class AnalyzerBondBondCorrelation : public AbstractAnalyzer
{

private:
    //! lower histrogram border
    const double low_dist_bound = -M_PI;
    //! upper histogram border
    const double upp_dist_bound = M_PI;
    //! number of bins for histogram
    const int n_bins = 30;
    //! length, as seperation distance, to calculate bbc for
    int correlation_length;
    //! typedef for the underlying container holding the monomers
    typedef typename IngredientsType::molecules_type molecules_type;
    //! reference to the complete system
    const IngredientsType& ingredients;
    //! bb is calculated for the groups in this vector
    std::vector<MonomerGroup<molecules_type> > groups;
    //! histogram that accummulates the bond bond angle values
    HistogramGeneralStatistik1D bbcorr;
    //! where to store analysis results
    std::string outputFile;
    //! calculate the bond bond distrubution  of the monomer group
    std::vector<double> cummulateBBC(const MonomerGroup<molecules_type>& group) const;
    //! analyze only after equilibration time has been reached
    uint32_t equilibrationTime;
protected:
    //! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
    void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

    //! constructor
    AnalyzerBondBondCorrelation(const IngredientsType& ing,
                                 std::string filename="AnalyzerBondBondCorrelation.dat",
                                 uint32_t equilibrationTime_=0,
                                 int correlation_length_ = 30.);
    //only used to make sure you initialize your groups before you do things
    bool initialized;
    //! destructor. does nothing
    virtual ~AnalyzerBondBondCorrelation(){}
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
AnalyzerBondBondCorrelation<IngredientsType>::AnalyzerBondBondCorrelation(const IngredientsType& ing, std::string filename, uint32_t equilibrationTime_, int correlation_length_)
: ingredients(ing), outputFile(filename), equilibrationTime(equilibrationTime_), correlation_length(correlation_length_)
{initialized=false;}

// initlaizer
template<class IngredientsType>
void AnalyzerBondBondCorrelation<IngredientsType>::initialize()
{

    // initialize histogram
    bbcorr = HistogramGeneralStatistik1D(low_dist_bound, upp_dist_bound, upp_dist_bound);
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
bool AnalyzerBondBondCorrelation<IngredientsType>::execute()
{

    //check if groups have been initialized. if not, exit and explain
    if(initialized==false)
    {
        std::stringstream errormessage;
        errormessage<<"AnalyzerBondBondCorrelation::execute()...analyzer not initialized\n"
        <<"Use AnalyzerBondBondCorrelation::initialize() or Taskmanager::init()\n";

        throw std::runtime_error(errormessage.str());
    }

    // if the equilibriation time has been reached, evalulate
    if(ingredients.getMolecules().getAge() >= equilibrationTime)
    {
        // loop through each monomer group (i.e. polymer)for(size_t n=0;n<groups.size();n++)
        for(size_t n=0;n<groups.size();n++) {
            // calculate the bond-bond angle for all pairs seperated by distance s
            for (int m = 1; m <= bbcorr.getNBins(); m ++) {
                std::vector<double> angles = cummulateBBC(groups[n]);
                // accumulate the angles into the histogram
            }

            // TODO implement cummulateBBC method
            // calculate cosTheta as a function of s, where s is the distance between i and j (abs |i-j|)
            std::vector<double> angles = cummulateBBC(groups[n]);
            for (int m = 0; m < angles.size(); m++) {
                bbcorr.addValue(angles[m], 1.);
            }
        }
    }

    return true;
}

template<class IngredientsType>
void AnalyzerBondBondCorrelation<IngredientsType>::cleanup()
{

    // open file stream
    std::ofstream file;
    file.open(outputFile.c_str());

    // TODO check if file is open, etc.

    // print results into a file
    double sum = 0.;
    for (int n = 0; n < bbcorr.getNBins(); n++) {
        sum += bbcorr.getFirstMomentInBin(n);
        file << bbcorr.getCenterOfBin(n) << ", " << bbcorr.getFirstMomentInBin(n) << std::endl;
    }
    file.close();
}

// cummulate bond bond angles for polymer chain
template<class IngredientsType>
std::vector<double> AnalyzerBondBondCorrelation<IngredientsType>::cummulateBBC(
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
            angle = M_PI * (bi*bj)/angle;
            // bbcorr.addValue(angle, 1.);
            angles[n] = angle;
        }
    return angles;
    }

#endif /*ANALYZER_BOND_BOND_CORRELATION_H*/


