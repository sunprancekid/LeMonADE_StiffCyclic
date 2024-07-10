#ifndef LEMONADE_ANALYZER_END_TO_END_DISTANCE_H
#define LEMONADE_ANALYZER_END_TO_END_DISTANCE_H

#include <string>
#include <iostream>
#include <cmath>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DistanceCalculation.h>

/*************************************************************************
 * definition of AnalyzerEndToEndDistance class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerEndToEndDistance
 *
 * @brief Analyzer for evaluating the squared radius of gyration Rg^2
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 *
 * @details This analyzer calculates the squared radius of gyration of either the
 * complete system (if no subgroup of monomers was specified), or the average of
 * the Rg^2 of a set of subgroups (if subgroups were specified using the provided
 * function void setMonomerGroups(std::vector<MonomerGroup<molecules_type> >) )
 * The analyzer calculates Rg^2 in every timestep and saves the time series
 * of this data to disk in the format
 * mcs <Re2e_x> <Re2_ey> <Re2e_z>  <Re2e_tot>
 * The default output filename is Re2eTimeSeries.dat and it can be changed by argument
 * to the constructor or by using the setter function provided.
 * If a more sophisticated grouping of monomers into groups is required, one can
 * also write a new analyzer, inheriting from AnalyzerEndToEndDistance, and
 * overwriting the initialize function.
 */

template < class IngredientsType > class AnalyzerEndToEndDistance : public AbstractAnalyzer
{

private:
    //! typedef for the underlying container holding the monomers
    typedef typename IngredientsType::molecules_type molecules_type;
    //! reference to the complete system
    const IngredientsType& ingredients;
    //! Rg2 is calculated for the groups in this vector
    std::vector<MonomerGroup<molecules_type> > groups;
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
    //! integer representing the furthest monomer used for end-to-end vector calculation
    //! boolean that determines if the minimum image convention should be applied to distance vector calculations
    // bool minImage;
protected:
    //! Set the groups to be analyzed. This function is meant to be used in initialize() of derived classes.
    void setMonomerGroups(std::vector<MonomerGroup<molecules_type> > groupVector){groups=groupVector;}
public:

    //! constructor
    AnalyzerEndToEndDistance(const IngredientsType& ing,
                             std::string filename="Re2eTimeSeries.dat",
                             uint32_t equilibrationTime_=0);

    //! destructor. does nothing
    virtual ~AnalyzerEndToEndDistance(){}
    //! Initializes data structures. Called by TaskManager::initialize()
    virtual void initialize();
    //! Calculates the Rg2 for the current timestep. Called by TaskManager::execute()
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
    //! monomers that are subjected to a force according to their attributes tag
    uint32_t indexMonomersForce[2];
    //! vector that total distance is projected onto
    VectorDouble3 projVec;
    //! boolean that determines if the distance should be calculated in terms of the projection vector
    bool calcProjVec;
    //! assigns vector that distance should be  projected onto
    void setE2EProjVec(double x, double y, double z) {projVec.setAllCoordinates(x, y, z); calcProjVec = true;};
    //! assigns projection vector as 3d double vector type
    void setE2EProjVec(VectorDouble3 v) {setE2EProjVec(v.getX(), v.getY(), v.getZ());};
    //! returns the projection vector
    VectorDouble3 getE2EProjVec() {return projVec;};

};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * @param filename output file name. defaults to "Re2eTimeSeries.dat".
 * */
template<class IngredientsType>
AnalyzerEndToEndDistance<IngredientsType>::AnalyzerEndToEndDistance(
    const IngredientsType& ing,
    std::string filename,
    uint32_t equilibrationTime_)
:ingredients(ing)
,Re2eTimeSeries(4,std::vector<double>(0))
,bufferSize(100)
,outputFile(filename)
,isFirstFileDump(true)
,equilibrationTime(equilibrationTime_)
,calcProjVec(false)
{
}

/**
 * @details fills all monomers into the groups vector as a single group,
 * if the group size is zero. Normally this should always be the case,
 * unless the groups were already set explicitly by using the setter function
 * setMonomerGroups
 * */
template< class IngredientsType >
void AnalyzerEndToEndDistance<IngredientsType>::initialize()
{
    // record the number of initial number of MCS steps
    startMCS=ingredients.getMolecules().getAge();
    //if no groups are set, use the complete system by default
    //groups can be set using the provided access function
    if(groups.size()==0){
        groups.push_back(MonomerGroup<molecules_type>(ingredients.getMolecules()));
        for(size_t n=0;n<ingredients.getMolecules().size();n++)
            groups[0].push_back(n);
    }
    //get the monomers with attributes subjected to force
    uint32_t n = 0;
    for (uint32_t i=0; i < ingredients.getMolecules().size();i++){
        if (ingredients.getMolecules()[i].getAttributeTag() != 0){
            // std::cout << ingredients.getMolecules()[i].getAttributeTag() << "\n";
            indexMonomersForce[n] = i;
            n += 1;
        }
    }
    // check if the constant force is on
    if (ingredients.isConstantForceOn()) {
        // if the constant force is on, set the projection vector as the force vector
        setE2EProjVec(ingredients.getForceVector());
        // std::cout << calcProjVec << std::endl;
        // exit(0);
    }
}

/**
 * @details Calculates the current Rg2, saves it in the
 * time series, and saves the time series to disk in regular intervals.
 * */
template< class IngredientsType >
bool AnalyzerEndToEndDistance<IngredientsType>::execute()
{
    VectorDouble3 Re2eComponents(0.0,0.0,0.0);

    // if the equilibriation time for the analyzer has not been reached yet
    if((ingredients.getMolecules().getAge()-startMCS) <= equilibrationTime)
        // skip the e2e vector calculation
        return true;

    for(size_t n=0;n<groups.size();n++)
    {
        //this vector will contain (Re2e_x, Re2e_y, Re2e_z)
        Re2eComponents+=calculateRe2eComponents(groups[n]);
        // cumulate components, distance vector
        Re2eTimeSeries[0].push_back(Re2eComponents.getX());
        Re2eTimeSeries[1].push_back(Re2eComponents.getY());
        Re2eTimeSeries[2].push_back(Re2eComponents.getZ());
        Re2eTimeSeries[3].push_back(std::sqrt(std::pow(Re2eComponents.getX(), 2)+std::pow(Re2eComponents.getY(),2)+std::pow(Re2eComponents.getZ(),2)));
        // std::cout << std::sqrt(std::pow(Re2eComponents.getX(), 2)+std::pow(Re2eComponents.getY(),2)+std::pow(Re2eComponents.getZ(),2)) << "\n";
        // exit(0);
    }

    MCSTimes.push_back(ingredients.getMolecules().getAge());
    //save to disk in regular intervals, only if equilibriation time has been reached
    if(MCSTimes.size()>=bufferSize)
        dumpTimeSeries();

    return true;
}


template<class IngredientsType>
void AnalyzerEndToEndDistance<IngredientsType>::cleanup()
{
    std::cout<<"AnalyzerEndToEndDistance::cleanup()...";
    //write the remaining data from the time series
    dumpTimeSeries();
    std::cout<<"done\n";
}

/**
 * @details Saves the current content of Re2eTimeSeries to the file
 * Re2eTimeSeries.dat. The output format is:
 * mcs Rg2_x Rg2_y Rg2_z Rg2_tot
 * */
template<class IngredientsType>
void AnalyzerEndToEndDistance<IngredientsType>::dumpTimeSeries()
{
    //fist make a single vector<vector<double> > for writing the results
    std::vector<std::vector<double> > resultsTimeseries=Re2eTimeSeries;
    resultsTimeseries.insert(resultsTimeseries.begin(),MCSTimes);

    //if it is written for the first time, include comment in the output file
    if(isFirstFileDump){
        std::stringstream commentTimeSeries;
        commentTimeSeries<<"Created by AnalyzerEndToEndDistance\n";
        commentTimeSeries<<"file contains time series of average Re2e over all analyzed groups\n";
        commentTimeSeries<<"format: mcs\t Re2e_X\t Re2e_Y\t Re2e_Z\t Re2e\n";

        ResultFormattingTools::writeResultFile(
            outputFile,
            ingredients,
            resultsTimeseries,
            commentTimeSeries.str()
        );

        isFirstFileDump=false;
    }
    //otherwise just append the new data
    else{
        ResultFormattingTools::appendToResultFile(outputFile,
                                                  resultsTimeseries);
    }
    //set all time series vectors back to zero size
    MCSTimes.resize(0);
    Re2eTimeSeries.resize(0);
    Re2eTimeSeries.resize(4,std::vector<double>(0));
}

/**
 * @details calculates the three components Rg^2_x, Rg^2_y,Rg^2_z and returns
 * them in a vector.
 * @return VectorDouble3 containing the components Rg^2_x, Rg^2_y,Rg^2_z, or (0.0,0.0,0.0) if group is empty)
 * @param group the monomer group of which the Rg2 is calculated
 * */
template<class IngredientsType>
VectorDouble3 AnalyzerEndToEndDistance<IngredientsType>::calculateRe2eComponents(
    const MonomerGroup<molecules_type>& group) const
    {
        // get the group size
        double group_size = group.size();

        //if group is empty, return zero vector
        if(group_size==0){
            return VectorDouble3(0.0,0.0,0.0);
        }

        // vector containing the difference between the first and last momomers in the group
        VectorDouble3 diff_vec;
        diff_vec.setX(group[indexMonomersForce[1]].getX() - group[indexMonomersForce[0]].getX());
        diff_vec.setY(group[indexMonomersForce[1]].getY() - group[indexMonomersForce[0]].getY());
        diff_vec.setZ(group[indexMonomersForce[1]].getZ() - group[indexMonomersForce[0]].getZ());
        if (calcProjVec) {
            double mag = (diff_vec * projVec) / (projVec * projVec);
            diff_vec = mag * projVec;
        }
        return diff_vec;
    }
    #endif
