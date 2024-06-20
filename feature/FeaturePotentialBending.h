/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2024 by
  o/.|.\o    E   nvironment    | Ron Dockhorn
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

#ifndef LEMONADE_FEATURE_FEATUREPOTENTIALBENDING_H
#define LEMONADE_FEATURE_FEATUREPOTENTIALBENDING_H

#include <tuple>

#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/io/AbstractRead.h>
#include <LeMonADE/io/AbstractWrite.h>
#include <LeMonADE/io/FileImport.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/updater/moves/MoveAddMonomerBase.h>
#include <LeMonADE/feature/FeatureBoltzmann.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>


/**
 * @class MonomerBendingBondInformation
 * @brief Extends monomers by a vector of indices of monomer forming consecutive bonds.
 * */

class MonomerBendingBondInformation
{
public:

	//! Standard constructor- initially the tag is set to NULL.
	MonomerBendingBondInformation(){}

	//! Getting the tag of the monomer.
	//TagType getAttributeTag() const {return tag;}
	
	std::vector<std::tuple<std::pair<size_t,size_t>, std::pair<size_t,size_t> > > getBendingBondInformation() const {return vBendingBonds;}

	/**
	 * @brief Setting the tag of the monomer with \para attr.
	 *
	 * @param attr
	 */
	//void setAttributeTag(TagType attribute){ tag=attribute;}
	
	void setBendingBondInformation(std::pair<size_t,size_t> firstBond, std::pair<size_t,size_t> secondBond)
    {
        vBendingBonds.push_back(std::make_tuple(firstBond, secondBond));
    }

private:
	 //! Private variable holding the tag. Default is NULL.
     //TagType tag;
     
     /** 
      * @brief Defines the angle for applying the bending potential\n
      * using the indices of the monomers. The first pair defines the first bond\n
      * and the second pair defines the second bond in the dot product for calculating cos(theta).\n
      * std::pair<int A,int B> => std::pair<index monomer A, index monomer B> == std::pair<initial point A of first bond, terminal point B of first bond>
      * std::pair<int C,int D> => std::pair<index monomer C, index monomer D> == std::pair<initial point C of second bond, terminal point D of second bond>
      * 
      *       D^
      *       /
      *      /  theta
      *    C/_________>
      *      A        B
      */ 
     std::vector<std::tuple<std::pair<size_t,size_t>, std::pair<size_t,size_t> > > vBendingBonds;
};


/*****************************************************************/
/**
 * @class ReadBendingBondInformation
 *
 * @brief Handles BFM-File-Reads \b #!bending_potential_bonds
 * @tparam IngredientsType Ingredients class storing all system information.
 */
template < class IngredientsType>
class ReadBendingBondInformation: public ReadToDestination<IngredientsType>
{
public:
  ReadBendingBondInformation(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
  virtual ~ReadBendingBondInformation(){}
  virtual void execute();
};


/*****************************************************************/
/**
 * @class WriteBendingBondInformation
 *
 * @brief Handles BFM-File-Write \b #!bending_potential_bonds
 * @tparam IngredientsType Ingredients class storing all system information.
 **/
template <class IngredientsType>
class WriteBendingBondInformation:public AbstractWrite<IngredientsType>
{
public:
	//! Only writes \b #!bending_potential_bonds into the header of the bfm-file.
  WriteBendingBondInformation(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
  virtual ~WriteBendingBondInformation(){}
  virtual void writeStream(std::ostream& strm);
};


/*****************************************************************/
/**
 * @class ReadBendingBondConstant
 *
 * @brief Handles BFM-File-Reads \b #!bending_potential_constant
 * @tparam IngredientsType Ingredients class storing all system information.
 **/

template < class IngredientsType>
class ReadBendingBondConstant: public ReadToDestination<IngredientsType>
{
public:
    ReadBendingBondConstant(IngredientsType& i):ReadToDestination<IngredientsType>(i){}
    virtual ~ReadBendingBondConstant(){}
    virtual void execute();
};

/*****************************************************************/
/**
 * @class WriteBendingBondConstant
 *
 * @brief Handles BFM-File-Write \b #!bending_potential_constant
 * @tparam IngredientsType Ingredients class storing all system information.
 **/

template < class IngredientsType>
class WriteBendingBondConstant:public AbstractWrite<IngredientsType>
{
public:
    WriteBendingBondConstant(const IngredientsType& i)
    :AbstractWrite<IngredientsType>(i){this->setHeaderOnly(true);}
        
    virtual ~WriteBendingBondConstant(){}
            	    	
    virtual void writeStream(std::ostream& strm);
};



/*****************************************************************/
/**
 * @class FeaturePotentialBending
 * @brief Extends vertex/monomer by the bending potential in the check move function.
 **/
class FeaturePotentialBending:public Feature
{
public:
    
  FeaturePotentialBending(): Bending_Potential_Constant(0.0){};
  virtual ~FeaturePotentialBending(){};
  
  //! This Feature requires a monomer_extensions.
  typedef LOKI_TYPELIST_1(MonomerBendingBondInformation) monomer_extensions;
  
  //the FeatureBoltzmann: adds a probability for the move -> done after all checks in the feature list 
  typedef LOKI_TYPELIST_1(FeatureBoltzmann) required_features_back;
  
  //FeatureExcludedVolumeSc needs to be in front, because the bonds have be behave well defined.
  typedef LOKI_TYPELIST_1(FeatureExcludedVolumeSc< >) required_features_front;
  

  //! Export the relevant functionality for reading bfm-files to the responsible reader object
  template<class IngredientsType>
  void exportRead(FileImport<IngredientsType>& fileReader);

  //! Export the relevant functionality for writing bfm-files to the responsible writer object
  template<class IngredientsType>
  void exportWrite(AnalyzerWriteBfmFile<IngredientsType>& fileWriter) const;

  //! For all unknown moves: this does nothing
  template<class IngredientsType> 
  bool checkMove(const IngredientsType& ingredients,const MoveBase& move) const{return true;};
    
  //! check for a MoveLocalSc 
  template<class IngredientsType> 
  bool checkMove(const IngredientsType& ingredients,MoveLocalSc& move) const;
        
  //! For all unknown moves: this does nothing
  template<class IngredientsType>
  void applyMove(IngredientsType& ing, const MoveBase& move){};

  //! Synchronize the system: Checks if bonds are valid.
  template<class IngredientsType>
  void synchronize(IngredientsType& ingredients);
  
  //! Set the bending potential constant
  void setBending_Potential_Constant(double value)
  {
      Bending_Potential_Constant = value;
  }
  
  //! Get the bending potential constant
  double getBending_Potential_Constant() const
  {
      return Bending_Potential_Constant;
  }
  
  
private:
    
  double Bending_Potential_Constant; //!< The constant as prefactor of the bending potential
   
};


/******************************************************************************
 * member implementations
 * ****************************************************************************/
template<class IngredientsType>
bool FeaturePotentialBending::checkMove(const IngredientsType& ingredients, MoveLocalSc& move) const
{
    //add the probability factor coming from this feature, then return true,
    //because the total probability is evaluated by FeatureBoltzmann at the end.
    
    const uint32_t monoIndex(move.getIndex());
    
    // we are calculating the bending potential VB=sum_{all angles} VB(theta_i)
    // and only the energetic difference dVB=VB(new config) - VB(old config)
    
    // get all bending bond information of the moving monomer
    // as excl. volume is present all bond vectors are well defined
    // for comparison of the bending potential have also a look on DOI: 10.3390/polym11020295
    auto  vBendingBonds = ingredients.getMolecules()[monoIndex].getBendingBondInformation();
    // is there any information about bonds in bending potential
    if(!vBendingBonds.empty())
    {
        double VB_old = 0.0;
        double VB_new = 0.0;
        
        // if so: loop over all
        for (auto bondinfo : vBendingBonds)
        {
            // get the position of the monomers in space
            VectorInt3 b1_offset = ingredients.getMolecules()[std::get<0>(bondinfo).first];
            VectorInt3 b1_terminal = ingredients.getMolecules()[std::get<0>(bondinfo).second];
            
            VectorInt3 b2_offset = ingredients.getMolecules()[std::get<1>(bondinfo).first];
            VectorInt3 b2_terminal = ingredients.getMolecules()[std::get<1>(bondinfo).second];
            
            // old bond
            VectorInt3 bondvector1_old = b1_terminal - b1_offset;
            VectorInt3 bondvector2_old = b2_terminal - b2_offset;;
            
            // Bending potential VB(theta_i)=KB*(1-cos(theta_i))=KB*[1-(a*b)/(|a|*|b|)]
            // see Martemyanova, J. A. / Stukan, M. R. / Ivanov, V. A. / Müller, M. / Paul, W. / Binder, K. 
            // Dense orientationally ordered states of a single semiflexible macromolecule: An expanded ensemble Monte Carlo simulation 
            // The Journal of Chemical Physics 2005, Vol. 122, No. 17 DOI:10.1063/1.1888525
            //
            // VB_old += Bending_Potential_Constant*(1.0-double(bondvector1_old*bondvector2_old)/(bondvector1_old.getLength()*bondvector2_old.getLength()));
            //
            
            // Bending potential VB(theta_i)=KB*(cos(theta_i)-1)²=KB*[(a*b)/(|a|*|b|)-1]²
            // see Ivanov, V. A. / Paul, W. / Binder, K. 
            // Finite chain length effects on the coil–globule transition of stiff-chain macromolecules: A Monte Carlo simulation 
            // The Journal of Chemical Physics 1998, Vol. 109, No. 13 DOI: 10.1063/1.477184
            //
            double tmp =(double(bondvector1_old*bondvector2_old)/(bondvector1_old.getLength()*bondvector2_old.getLength())-1.0);
            VB_old += Bending_Potential_Constant*(tmp*tmp);
            //
            
            // Bending potential VB(theta_i)=-KB*cos(theta_i - theta0)=>(theta0=0):VB(theta_i)=-KB*[(a*b)/(|a|*|b|)]
            // Bates, Martin A. 
            // Coarse grained models for flexible liquid crystals: Parameterization of the bond fluctuation model 
            // The Journal of Chemical Physics 2004, Vol. 120, No. 4 DOI: 10.1063/1.1634551
            //
            //VB_old += -Bending_Potential_Constant*(double(bondvector1_old*bondvector2_old)/(bondvector1_old.getLength()*bondvector2_old.getLength()));
            //
            
            // now calculate the attempted new energetic contriubution
            // check for the moving monomer
            if(std::get<0>(bondinfo).first == monoIndex)
                b1_offset = ingredients.getMolecules()[monoIndex] + move.getDir();
            
            if(std::get<0>(bondinfo).second == monoIndex)
                b1_terminal = ingredients.getMolecules()[monoIndex] + move.getDir();
            
            if(std::get<1>(bondinfo).first == monoIndex)
                b2_offset = ingredients.getMolecules()[monoIndex] + move.getDir();
            
            if(std::get<1>(bondinfo).second == monoIndex)
                b2_terminal = ingredients.getMolecules()[monoIndex] + move.getDir();
            
            // new bond
            VectorInt3 bondvector1_new = b1_terminal - b1_offset;
            VectorInt3 bondvector2_new = b2_terminal - b2_offset;
            
            // Bending potential VB(theta_i)=KB*(1-cos(theta_i))=KB*[1-(a*b)/(|a|*|b|)]
            // see Martemyanova, J. A. / Stukan, M. R. / Ivanov, V. A. / Müller, M. / Paul, W. / Binder, K. 
            // Dense orientationally ordered states of a single semiflexible macromolecule: An expanded ensemble Monte Carlo simulation 
            // The Journal of Chemical Physics 2005, Vol. 122, No. 17 DOI:10.1063/1.1888525
            //
            // VB_new += Bending_Potential_Constant*(1.0-double(bondvector1_new*bondvector2_new)/(bondvector1_new.getLength()*bondvector2_new.getLength()));
            //
            
            // Bending potential VB(theta_i)=KB*(cos(theta_i)-1)²=KB*[(a*b)/(|a|*|b|)-1]²
            // see Ivanov, V. A. / Paul, W. / Binder, K. 
            // Finite chain length effects on the coil–globule transition of stiff-chain macromolecules: A Monte Carlo simulation 
            // The Journal of Chemical Physics 1998, Vol. 109, No. 13 DOI: 10.1063/1.477184
            //
            double tmp_b = (double(bondvector1_new*bondvector2_new)/(bondvector1_new.getLength()*bondvector2_new.getLength())-1.0);
            VB_new += Bending_Potential_Constant*(tmp_b*tmp_b);
            //
            
            // Bending potential VB(theta_i)=-KB*cos(theta_i - theta0)=>(theta0=0):VB(theta_i)=-KB*[(a*b)/(|a|*|b|)]
            // Bates, Martin A. 
            // Coarse grained models for flexible liquid crystals: Parameterization of the bond fluctuation model 
            // The Journal of Chemical Physics 2004, Vol. 120, No. 4 DOI: 10.1063/1.1634551
            //
            // VB_new += -Bending_Potential_Constant*(double(bondvector1_new*bondvector2_new)/(bondvector1_new.getLength()*bondvector2_new.getLength()));
            //
            
        }
       
       //std::cout << VB_new << " - " << VB_old << " = " << (VB_new-VB_old) << std::endl;
       // apply the probability to move
       double prob=exp(-(VB_new-VB_old));
       move.multiplyProbability(prob);
    }
    
    return true;
}


/******************************************************************************/
/**
 * @fn void FeatureConnectionSc ::synchronize(IngredientsType& ingredients)
 * @brief Synchronizes the lattice occupation with the rest of the system
 * by calling the private function fillLattice.
 *
 * @param ingredients A reference to the IngredientsType - mainly the system.
 */
/******************************************************************************/
template<class IngredientsType>
void FeaturePotentialBending::synchronize(IngredientsType& ingredients)
{
    
    std::cout << "FeaturePotentialBending::Checking bond information for bending...\n";
    // all bonds for bending calculation should also be connected bonds
   
    for(size_t i = 0 ; i < ingredients.getMolecules().size(); i++ ){
        
        // is there any information about bonds in bending potential
        if(!ingredients.getMolecules()[i].getBendingBondInformation().empty())
        {
            // if so: loop over all
            for (auto bondinfo : ingredients.getMolecules()[i].getBendingBondInformation()){
                
                // uncomment if you need more output
                // std::cout << i << ":" << std::get<0>(bondinfo).first << ", " << std::get<0>(bondinfo).second << ", " << std::get<1>(bondinfo).first << ", " << std::get<1>(bondinfo).second << '\n';
                
                // check if they are connected
                if ( !ingredients.getMolecules().areConnected(std::get<0>(bondinfo).first, std::get<0>(bondinfo).second) )
                {
                    throw std::runtime_error("********** FeaturePotentialBending::synchronize(): bending bond angle calculation of unconnected bonds ******************");
                }
                
                if ( !ingredients.getMolecules().areConnected(std::get<1>(bondinfo).first, std::get<1>(bondinfo).second) )
                {
                    throw std::runtime_error("********** FeaturePotentialBending::synchronize(): bending bond angle calculation of unconnected bonds ******************");
                }
            }
        }
        
    }
    std::cout << "done without error" << std::endl;
}


/**
 * @details The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type FileImport. The export of the Reads is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Read-In Commands:
 * * !attributes
 *
 * @param fileReader File importer for the bfm-file
 * @param destination List of Feature to write-in from the read values.
 * @tparam IngredientsType Features used in the system. See Ingredients.
 **/
template<class IngredientsType>
void FeaturePotentialBending::exportRead(FileImport< IngredientsType >& fileReader)
{
  fileReader.registerRead("#!bending_potential_bonds",new ReadBendingBondInformation<IngredientsType>(fileReader.getDestination()));
  fileReader.registerRead("#!bending_potential_constant",new ReadBendingBondConstant<IngredientsType>(fileReader.getDestination()));
}


/**
 * The function is called by the Ingredients class when an object of type Ingredients
 * is associated with an object of type AnalyzerWriteBfmFile. The export of the Writes is thus
 * taken care automatically when it becomes necessary.\n
 * Registered Write-Out Commands:
 * * !attributes
 *
 * @param fileWriter File writer for the bfm-file.
 */
template<class IngredientsType>
void FeaturePotentialBending::exportWrite(AnalyzerWriteBfmFile< IngredientsType >& fileWriter) const
{
  fileWriter.registerWrite("#!bending_potential_bonds",new WriteBendingBondInformation<IngredientsType>(fileWriter.getIngredients_()));
  fileWriter.registerWrite("#!bending_potential_constant", new WriteBendingBondConstant<IngredientsType>(fileWriter.getIngredients_()));
}



/**
 * @brief Executes the reading routine to extract \b #!bending_potential_bonds 
 *
 * @throw <std::runtime_error> attributes and identifier could not be read.
 **/
template < class IngredientsType >
void ReadBendingBondInformation<IngredientsType>::execute()
{
  //some variables used during reading
  //counts the number of attribute lines in the file
  int nAttributes=0;
  
  // file layout for linear chain with 5 monomers (5-2=3 angels):
  // #!bending_potential_bonds
  // 1:1-2;2-3
  // 2:1-2;2-3
  // 2:2-3;3-4
  // 3:1-2;2-3
  // 3:2-3;3-4
  // 3:3-4;4-5
  // 4:2-3;3-4
  // 4:3-4;4-5
  // 5:3-4;4-5
  int startIndex;
  int startIndexFirstBond, endIndexFirstBond;
  int startIndexSecondBond, endIndexSecondBond;
  
  //TagType attribute;
  //contains the latest line read from file
  std::string line;
  //used to reset the position of the get pointer after processing the command
  std::streampos previous;
  //for convenience: get the input stream
  std::istream& source=this->getInputStream();
  //for convenience: get the set of monomers
  typename IngredientsType::molecules_type& molecules=this->getDestination().modifyMolecules();

  //go to next line and save the position of the get pointer into streampos previous
  getline(source,line);
  previous=(source).tellg();

  //read and process the lines containing the bond vector definition
  getline(source,line);

  while(!line.empty() && !((source).fail())){

    //stop at next Read and set the get-pointer to the position before the Read
    if(this->detectRead(line)){
      (source).seekg(previous);
      break;
    }

    //initialize stringstream with content for ease of processing
    std::stringstream stream(line);

    //read monomer index
    stream>>startIndex;

    //throw exception, if extraction fails
    if(stream.fail()){
      std::stringstream messagestream;
      messagestream<<"ReadBendingBondInformation<IngredientsType>::execute()\n"
		   <<"Could not read first index in attributes line "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt "-"
    if(!this->findSeparator(stream,':')){

    	std::stringstream messagestream;
      messagestream<<"ReadBendingBondInformation<IngredientsType>::execute()\n"
		   <<"Wrong definition of attributes\nCould not find separator \":\" "
		   <<"in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }

    //read first bond first monomer index == offset first bond, throw exception if extraction fails
    stream>>startIndexFirstBond;
    
    //throw exception, if extraction fails
    if(stream.fail()){
    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Could not read second index in attributes line "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt "-"
    if(!this->findSeparator(stream,'-')){

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Wrong definition of attributes\nCould not find separator \"-\" "
		   <<"in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }
    
    //read first bond second monomer index == terminal first bond, throw exception if extraction fails
    stream>>endIndexFirstBond;
    
    //throw exception, if extraction fails
    if(stream.fail()){
    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Could not read second index in attributes line "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());
    }
    
    //throw exception, if next character isnt ";"
    if(!this->findSeparator(stream,';')){

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Wrong definition of attributes\nCould not find separator \";\" "
		   <<"in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }
    
    //read second bond first monomer index == offset second bond, throw exception if extraction fails
    stream>>startIndexSecondBond;
    
    //throw exception, if extraction fails
    if(stream.fail()){
    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Could not read second index in attributes line "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());
    }

    //throw exception, if next character isnt "-"
    if(!this->findSeparator(stream,'-')){

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Wrong definition of attributes\nCould not find separator \"-\" "
		   <<"in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }
    
    //read second bond second monomer index == terminal second bond, throw exception if extraction fails
    stream>>endIndexSecondBond;
    
    //throw exception, if extraction fails
    if(stream.fail()){
    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"Could not read second index in attributes line "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());
    }
    
    //if extraction worked, save the information
    if(!stream.fail()){

      //save attributes
        //use n-1 as index, because bfm-files start counting indices at 1 (not 0)
        std::pair<size_t, size_t> firstBond{startIndexFirstBond-1, endIndexFirstBond-1};
        std::pair<size_t, size_t> secondBond{startIndexSecondBond-1, endIndexSecondBond-1};
        molecules[startIndex-1].setBendingBondInformation(firstBond, secondBond);
        
      nAttributes++;
      getline((source),line);

    }
    //otherwise throw an exception
    else{

    	std::stringstream messagestream;
      messagestream<<"ReadAttributes<IngredientsType>::execute()\n"
		   <<"could not read attribute in attribute definition no "<<nAttributes+1;
      throw std::runtime_error(messagestream.str());

    }
  }
}


//! Executes the routine to write \b #!bending_potential_bonds.
template < class IngredientsType>
void WriteBendingBondInformation<IngredientsType>::writeStream(std::ostream& strm)
{
  //for all output the indices are increased by one, because the file-format
  //starts counting indices at 1 (not 0)
    
  // file layout for linear chain with 5 monomers (5-2=3 angels):
  // #!bending_potential_bonds
  // 1:1-2;2-3
  // 2:1-2;2-3
  // 2:2-3;3-4
  // 3:1-2;2-3
  // 3:2-3;3-4
  // 3:3-4;4-5
  // 4:2-3;3-4
  // 4:3-4;4-5
  // 5:3-4;4-5

  //write bfm command
  strm<<"#!bending_potential_bonds\n";
  //get reference to monomers
  const typename IngredientsType::molecules_type& molecules=this->getSource().getMolecules();

  size_t nMonomers=molecules.size();
  //attribute blocks begin with startIndex
  size_t startIndex=0;
  //counter varable
  size_t n=0;
  //attribute to be written (updated in loop below)
  //TagType attribute=molecules[0].getAttributeTag();

  //write bending bond information for each monomer (blockwise)
  while(n<nMonomers){
      // is there any information
    if(!molecules[n].getBendingBondInformation().empty())
    {
        for (auto i : molecules[n].getBendingBondInformation()){
            // don't forget: indices in file are increased by one
            strm << n+1 << ":" <<(std::get<0>(i).first+1)<<"-"<<(std::get<0>(i).second+1)<< ";" <<(std::get<1>(i).first+1)<<"-"<<(std::get<1>(i).second+1) <<std::endl;
        }
    }
    n++;
  }
  //write final line
  strm<<std::endl;

}



template < class IngredientsType >
void ReadBendingBondConstant<IngredientsType>::execute()
{
    std::cout<<"reading BendingBondConstant...";
    
    double amplitudeForce = 0.0;
    IngredientsType& ingredients=this->getDestination();
    std::istream& source=this->getInputStream();
    
    std::string line;
    getline(source,line);
    amplitudeForce = atof(line.c_str());
    std::cout << "#!bending_potential_constant=" << (amplitudeForce) << std::endl;
    
    ingredients.setBending_Potential_Constant(amplitudeForce);
}



template < class IngredientsType>
void WriteBendingBondConstant<IngredientsType>::writeStream(std::ostream& stream)
{
    stream<<"#!bending_potential_constant=" << (this->getSource().getBending_Potential_Constant()) << std::endl<< std::endl;
}


#endif /* LEMONADE_FEATURE_FEATUREPOTENTIALBENDING_H */
