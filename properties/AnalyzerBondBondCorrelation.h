
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
    //! lower histrogram border for bond bond angles
    const double low_bond_angle_bound = -1.;
    //! upper histogram border for bond bond angles
    const double upp_bond_angle_bound = 1.;
    //! lower histogram border for bond length
    const double low_bond_len_bound = 0.;
    //! upper histogram border for bond length
    const double upp_bond_len_bound = 4.;
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
    std::vector<HistogramGeneralStatistik1D> bbcorr;
    //! where to store analysis results
    std::string outputFile;
    //! calculate the bond bond distrubution  of the monomer group
    std::vector<double> cummulateBBC(const MonomerGroup<molecules_type>& group, int s) const;
    //! analyze only after equilibration time has been reached
    uint32_t equilibrationTime;
    //! histogram bin tolerance
    double TOL = 0.01;
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
    bbcorr.resize(correlation_length + 1); // resize the array to the correlation length + 1
    for (int n = 0; n < bbcorr.size(); n++) {
        if (n == 0) {
            // for the first histogram, initialize between the min and max bond length
            bbcorr[n] = HistogramGeneralStatistik1D(low_bond_len_bound  - TOL, upp_bond_len_bound + TOL, n_bins);
        }
        else {
            // initialize each element in the array with a new histogram
            bbcorr[n] = HistogramGeneralStatistik1D(low_bond_angle_bound  - TOL, upp_bond_angle_bound + TOL, n_bins);
        }
    }
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
        // loop through each monomer group (i.e. polymer)
        for(size_t n = 0; n < groups.size(); n++) {
            // calculate the bond-bond angle for all pairs seperated by distance s
            for (int m = 0; m < bbcorr.size(); m ++) {
                // std::cout << " m = " << m << std::endl;
                std::vector<double> angles = cummulateBBC(groups[n], m);
                // accumulate the angles into the histogram corresponding to the correlation length
                for (int k = 0; k < angles.size(); k++) {
                    // std::cout << " k = " << k << std::endl;
                    bbcorr[m].addValue(angles[k], 1.);
                }
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
    for (int n = 0; n <= correlation_length; n++) {
        file << (n) << ", " << bbcorr[n].getHistAverage() << std::endl;
    }
    file.close();
}

// cummulate bond bond angles for polymer chain
template<class IngredientsType>
std::vector<double> AnalyzerBondBondCorrelation<IngredientsType>::cummulateBBC(
    const MonomerGroup<molecules_type>& group, int s) const
    {
        // get the group size
        double group_size = group.size();

        // list of angles
        std::vector<double> angles;

        // if the group size is less than the correlation length + 1
        if(group_size<=s+1){
            // return an empty vector as it is there are no bonds seperated by this distance
            return angles;
        }
        if (s == 0) {
            // for the zeroth correlation length
            // accumulate the length of each bond vector
            angles.resize(group_size - 1);
            for (int n = 0; n < group_size - 1; n++) {
                // get the vector corresponding to the nth bond
                VectorDouble3 bi;
                bi.setX(group[n+1].getX() - group[n].getX());
                bi.setY(group[n+1].getY() - group[n].getY());
                bi.setZ(group[n+1].getZ() - group[n].getZ());
                // calculate the length
                double len = bi*bi;
                len = sqrt(len);
                angles[n] = len;
            }
        } else {
            // accumulate the angle between bonds seperated by distance s
            // the number of unique bond pairs seperated by seperation distance s is
            int N_s = group_size - s - 1;
            angles.resize(N_s);
            // loop through each unqiue bond pairs corresponding to the sepereration distance
            for (int n = 0; n < N_s; n++){
                // get the vectors representing the first and second bonds
                VectorDouble3 bi;
                bi.setX(group[n+1].getX() - group[n].getX());
                bi.setY(group[n+1].getY() - group[n].getY());
                bi.setZ(group[n+1].getZ() - group[n].getZ());
                VectorDouble3 bj;
                bj.setX(group[n+s+1].getX() - group[n+s].getX());
                bj.setY(group[n+s+1].getY() - group[n+s].getY());
                bj.setZ(group[n+s+1].getZ() - group[n+s].getZ());
                // calculate the angle between the two bonds
                double angle;
                angle = bi*bi;
                angle *= bj*bj;
                angle = sqrt(angle);
                angle = (bi*bj)/angle;
                angles[n] = angle;
            }
        }
    return angles;
    }

#endif /*ANALYZER_BOND_BOND_CORRELATION_H*/


