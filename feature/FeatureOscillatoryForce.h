/*--------------------------------------------------------------------------------
 o oo      L   attice-based  |*
 o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
 oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
 o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
 ooo                        |
 ----------------------------------------------------------------------------------

 This file is part of LeMonADE.

 LeMonADE is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 LeMonADE is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

 --------------------------------------------------------------------------------*/

#ifndef LEMONADE_OSCILLATORYFORCE_H
#define LEMONADE_OSCILLATORYFORCE_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveLocalSc.h>
#include <LeMonADE/updater/moves/MoveLocalScDiag.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/io/FileImport.h>
#include <bits/stdc++.h> // used to parse comma seperated force vector

/*****************************************************************************/
/**
 * @file
 * @date   2021/03/18
 * @author Toni & Bodonki
 *
 * @class FeatureOscillatoryForce
 * @brief This Feature applies a force on certain monomers.
 * @todo Substitute the FeatureAttributes<> by appropriate monomer extension
 * */
/*****************************************************************************/




class FeatureOscillatoryForce:public Feature
{
public:

    // constructor: assign default force vector and base force
    FeatureOscillatoryForce(): ConstantForceOn(false){setForceVector(1., 0., 0.); setBaseForce(0.);};
    virtual ~FeatureOscillatoryForce(){};

    //the FeatureBoltzmann: adds a probability for the move
    //the FeatureAttributes<>: labels the monomers for beeing force sensitive
    typedef LOKI_TYPELIST_2(FeatureBoltzmann,FeatureAttributes<>) required_features_back;

    //! For all unknown moves: this does nothing
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const{return true;};

    //! check for a MoveLocalSc
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const;

    //! check for a MoveLocalScDiag
    template<class IngredientsType>
    bool checkMove(const IngredientsType& ingredients,MoveLocalScDiag& move) const;

    //! force vector setter
    void setForceVector(double x, double y, double z) {
        forceVec.setAllCoordinates(x, y, z);
        // assign the normalization factor
        normfact = sqrt(forceVec * forceVec);
        // if a force has already been assigned, recalculate the probabilities
        setBaseForce(Base_Force);
    }

    //! return vector representing force magnitude and direction
    VectorDouble3 getForceVector() const {return forceVec;}

    //! set the strength of the force
    void setBaseForce(double baseForce){
        Base_Force = baseForce;
        prob.setAllCoordinates(exp(-(Base_Force / normfact) * forceVec.getX()), exp(-(Base_Force / normfact) * forceVec.getY()), exp(-(Base_Force / normfact) * forceVec.getZ()));
    }

    //! set force on or off
    void setForceOn(bool forceOn){
        // set boolean
        ConstantForceOn = forceOn;
    }

    //! get strength of the force
    double getBaseForce() const {return Base_Force;}

    //! applies a force on the monomers or not
    bool isConstantForceOn() const {return ConstantForceOn;}

    //! Export the relevant functionality for reading bfm-files to the responsible reader object
    template <class IngredientsType>
    void exportRead(FileImport <IngredientsType>& fileReader);

    //! Export the relevant functionality for writing bfm-files to the responsible writer object
    template <class IngredientsType>
    void exportWrite(AnalyzerWriteBfmFile <IngredientsType>& fileWriter) const;
private:

    //! vector representing force orientation
    VectorDouble3 forceVec;
    //! Force is On (True) or Off (False)
    bool ConstantForceOn;
    //! magnitude of base force (without oscilation / amplitude)
    double Base_Force;
    //! boolean determining if force oscilation is on or off
    // TODO add oscilatory force (force base, amplitude and period are specified)
    //! double representing oscilatory force period (MCS)
    //! double represnting force amplitude
    // TODO add method for turning on osciallatory force, setting period and amplitude
    // TODO add method for returning osciallatory force period
    // TODO add method for returning osciallatory force amplitude
    // TODO add method for checking if osciallatory force is on or off
    //!probability
    //! normalization factor for calculating acceptance probability in each dimension
    double normfact;
    //! move probability in each simulation dimension according to force magnitude and vector
    VectorDouble3 prob;

};

////////////////////////////////////////////////////////////////////////////////
//////////define member functions //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<class IngredientsType>
bool FeatureOscillatoryForce::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) const
{
    if(ConstantForceOn){

        const uint32_t monoIndex(move.getIndex());
        const int32_t tag(ingredients.getMolecules()[monoIndex].getAttributeTag());
        const int32_t dx(move.getDir().getX());
        const int32_t dy(move.getDir().getY());
        const int32_t dz(move.getDir().getZ());
        //Metropolis: zeta = exp (-dV)
        //dV=f*dr
        // positive force applied on attribute 4 in positive forceVec direction
        // negative force applied on attribute 5 in negative forceVec direction
        // NOTE :: algorithm assumes force vector is positive for all cartessian coordinates (this would be a problem for a rotating force vector)

        if ( dx == 1 ) {
            // positive x-direction
            if ( tag == 4 ) {
                // POSITIVE charge is MORE likely to move in the SAME direction as the force vector
                move.multiplyProbability(prob.getX());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is LESS likely to move in the SAME direction as the force vector
                move.multiplyProbability(1./prob.getX());
                return true;
            }
        } else if ( dx == -1 ) { // negative x-direction
            if ( tag == 4 ) {
                // POSITIVE charge is LESS likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(1./prob.getX());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is MORE likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(prob.getX());
                return true;
            }
        } else if ( dy == 1 ) {
            // positive y-direction
            if ( tag == 4 ) {
                // POSITIVE charge is MORE likely to move in the SAME direction as the force vector
                move.multiplyProbability(prob.getY());
                return true;
            } else if (tag == 5 ) {
                // NEGATIVE charge is LESS likely to move in the SAME direction as the force vector
                move.multiplyProbability(1./prob.getY());
                return true;
            }
        } else if ( dy == -1 ) {
            // negative y-direction
            if ( tag == 4 ) {
                // POSITIVE charge is LESS likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(1./prob.getY());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is MORE likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(prob.getY());
                return true;
            }
        } else if ( dz == 1 ) {
            // positive z-direction
            if ( tag == 4 ) {
                // POSITIVE charge is MORE likely to move in the SAME direction as the force vector
                move.multiplyProbability(prob.getZ());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is LESS likely to move in the SAME direction as the force vector
                move.multiplyProbability(1./prob.getZ());
                return true;
            }
        } else if ( dz == -1 ) {
            // negative z-direction
            if ( tag == 4 ) {
                // POSITIVE charge is LESS likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(1./prob.getZ());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is MORE likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(prob.getZ());
                return true;
            }
        }
    }
    return true;
}

template<class IngredientsType>
bool FeatureOscillatoryForce::checkMove(const IngredientsType& ingredients, MoveLocalScDiag& move) const
{
    if(ConstantForceOn){
        const uint32_t monoIndex(move.getIndex());
        const int32_t tag(ingredients.getMolecules()[monoIndex].getAttributeTag());
        const int32_t dx(move.getDir().getX());
        const int32_t dy(move.getDir().getY());
        const int32_t dz(move.getDir().getZ());
        //Metropolis: zeta = exp (-dV)
        //dV=f*dr
        // positive force applied on attribute 4 in positive forceVec direction
        // negative force applied on attribute 5 in negative forceVec direction
        // NOTE :: algorithm assumes force vector is positive for all cartessian coordinates (this would be a problem for a rotating force vector)

        if ( dx == 1 ) {
            // positive x-direction
            if ( tag == 4 ) {
                // POSITIVE charge is MORE likely to move in the SAME direction as the force vector
                move.multiplyProbability(prob.getX());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is LESS likely to move in the SAME direction as the force vector
                move.multiplyProbability(1./prob.getX());
                return true;
            }
        } else if ( dx == -1 ) { // negative x-direction
            if ( tag == 4 ) {
                // POSITIVE charge is LESS likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(1./prob.getX());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is MORE likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(prob.getX());
                return true;
            }
        } else if ( dy == 1 ) {
            // positive y-direction
            if ( tag == 4 ) {
                // POSITIVE charge is MORE likely to move in the SAME direction as the force vector
                move.multiplyProbability(prob.getY());
                return true;
            } else if (tag == 5 ) {
                // NEGATIVE charge is LESS likely to move in the SAME direction as the force vector
                move.multiplyProbability(1./prob.getY());
                return true;
            }
        } else if ( dy == -1 ) {
            // negative y-direction
            if ( tag == 4 ) {
                // POSITIVE charge is LESS likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(1./prob.getY());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is MORE likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(prob.getY());
                return true;
            }
        } else if ( dz == 1 ) {
            // positive z-direction
            if ( tag == 4 ) {
                // POSITIVE charge is MORE likely to move in the SAME direction as the force vector
                move.multiplyProbability(prob.getZ());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is LESS likely to move in the SAME direction as the force vector
                move.multiplyProbability(1./prob.getZ());
                return true;
            }
        } else if ( dz == -1 ) {
            // negative z-direction
            if ( tag == 4 ) {
                // POSITIVE charge is LESS likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(1./prob.getZ());
                return true;
            } else if ( tag == 5 ) {
                // NEGATIVE charge is MORE likely to move in the OPPOSITE direction of the force vector
                move.multiplyProbability(prob.getZ());
                return true;
            }
        }
    }
    return true;
}


/*****************************************************************/
/**
 * @class ReadForceFieldVector
 *
 * @brief Handles BFM-File-Reads \b #!force_field_vector
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadForceFieldVector: public ReadToDestination<IngredientsType>
{
public:
    ReadForceFieldVector(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadForceFieldVector(){}
    virtual void execute();
};

template<class IngredientsType>
void ReadForceFieldVector<IngredientsType>::execute()
{
    std::cout<<"reading ForceFieldVector...";

    // bool fieldIsOn = false;
    IngredientsType& ingredients=this->getDestination();
    std::istream& source=this->getInputStream();

    std::string line; // line parsed from config file
    vector<string> v; // contains each vector component
    getline(source,line); // parse line from file
    std::stringstream ss(line);
    while (ss.good()) {
        // loop through string
        string substr;
        getline(ss, substr, ',');
        v.push_back(substr);
    }
    // set force vec, report to user
    double x = std::stod(v[0]);
    double y = std::stod(v[1]);
    double z = std::stod(v[2]);
    std::cout << "#!force_field_vector=" << x << "," << y << "," << z << std::endl;

    ingredients.setForceVector(x, y, z);
}

/*****************************************************************/
/**
 * @class WriteForceFieldVector
 *
 * @brief Handles BFM-File-Write \b #!force_field_vector
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteForceFieldVector:public AbstractWrite<IngredientsType>
{
public:
    WriteForceFieldVector(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(false);}

    virtual ~WriteForceFieldVector(){}

    virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteForceFieldVector<IngredientsType>::writeStream(std::ostream& stream)
{
    VectorDouble3 fv = this->getSource().getForceVector();
    stream<<"#!force_vector=" << fv.getX() << "," << fv.getY() << "," << fv.getZ() << std::endl<< std::endl;
}

/*****************************************************************/
/**
 * @class ReadForceFieldOn
 *
 * @brief Handles BFM-File-Reads \b #!force_field_on
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadForceFieldOn: public ReadToDestination<IngredientsType>
{
public:
    ReadForceFieldOn(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadForceFieldOn(){}
    virtual void execute();
};

template<class IngredientsType>
void ReadForceFieldOn<IngredientsType>::execute()
{
    std::cout<<"reading FieldIsOn...";

    bool fieldIsOn = false;
    IngredientsType& ingredients=this->getDestination();
    std::istream& source=this->getInputStream();

    std::string line;
    getline(source,line);
    fieldIsOn = atoi(line.c_str());
    std::cout << "#!force_field_on=" << (fieldIsOn? "True" : " False" ) << std::endl;

    ingredients.setForceOn(fieldIsOn);
}

/*****************************************************************/
/**
 * @class WriteForceFieldOn
 *
 * @brief Handles BFM-File-Write \b #!force_field_on
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteForceFieldOn:public AbstractWrite<IngredientsType>
{
public:
    WriteForceFieldOn(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(false);}

    virtual ~WriteForceFieldOn(){}

    virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteForceFieldOn<IngredientsType>::writeStream(std::ostream& stream)
{
    stream<<"#!force_field_on=" << (this->getSource().isConstantForceOn() ? "1" : "0") << std::endl<< std::endl;
}


/*****************************************************************/
/**
 * @class ReadAmplitudeForce
 *
 * @brief Handles BFM-File-Reads \b #!force_amplitude
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template < class IngredientsType>
class ReadAmplitudeForce: public ReadToDestination<IngredientsType>
{
public:
    ReadAmplitudeForce(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadAmplitudeForce(){}
    virtual void execute();
};

template<class IngredientsType>
void ReadAmplitudeForce<IngredientsType>::execute()
{
    std::cout<<"reading AmplitudeForce...";

    double amplitudeForce = 0.0;
    IngredientsType& ingredients=this->getDestination();
    std::istream& source=this->getInputStream();

    std::string line;
    getline(source,line);
    amplitudeForce = atof(line.c_str());
    std::cout << "#!force_amplitude=" << (amplitudeForce) << std::endl;

    ingredients.setBaseForce(amplitudeForce);
}


/*****************************************************************/
/**
 * @class WriteAmplitudeForce
 *
 * @brief Handles BFM-File-Write \b #!force_amplitude
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteAmplitudeForce:public AbstractWrite<IngredientsType>
{
public:
    WriteAmplitudeForce(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(false);}

    virtual ~WriteAmplitudeForce(){}

    virtual void writeStream(std::ostream& strm);
};

template<class IngredientsType>
void WriteAmplitudeForce<IngredientsType>::writeStream(std::ostream& stream)
{
    stream<<"#!force_amplitude=" << (this->getSource().getBaseForce()) << std::endl<< std::endl;
}



/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * #!shear_field_on
 * * #!force_amplitude
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template<class IngredientsType>
void FeatureOscillatoryForce::exportRead(FileImport< IngredientsType >& fileReader)
{
    fileReader.registerRead("#!force_field_on", new ReadForceFieldOn<FeatureOscillatoryForce>(*this));
    fileReader.registerRead("#!force_amplitude", new ReadAmplitudeForce<FeatureOscillatoryForce>(*this));
    fileReader.registerRead("#!force_vector", new ReadForceFieldVector<FeatureOscillatoryForce>(*this));

}

/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * #!shear_field_on
 * * #!force_amplitude
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeatureOscillatoryForce::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
    fileWriter.registerWrite("#!force_field_on",new WriteForceFieldOn<FeatureOscillatoryForce>(*this));
    fileWriter.registerWrite("#!force_amplitude", new WriteAmplitudeForce<FeatureOscillatoryForce>(*this));
    fileWriter.registerWrite("#!force_vector", new WriteForceFieldVector<FeatureOscillatoryForce>(*this));
}




#endif /*LEMONADE_OSCILLATORYFORCE_H*/
