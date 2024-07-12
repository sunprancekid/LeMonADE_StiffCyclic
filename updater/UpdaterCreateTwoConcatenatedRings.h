/*--------------------------------------------------------------------------------
    ooo      L   attice-based  | LeMonADE: An Open Source Implementation of the
  o\.|./o    e   xtensible     |           Bond-Fluctuation-Model for Polymers
 o\.\|/./o   Mon te-Carlo      |           
oo---0---oo  A   lgorithm and  | StiffCyclic-simulations
 o/./|\.\o   D   evelopment    | Copyright (C) 2024 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
    ooo      StiffCyclic       | Ron Dockhorn
----------------------------------------------------------------------------------

This file is part of LeMonADE and LeMonADE_StiffCyclic extension.

LeMonADE and LeMonADE_StiffCyclic extension is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE  and LeMonADE_StiffCyclic extension is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE and LeMonADE_StiffCyclic extension. If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_UPDATER_CREATE_TWO_CONCATENATED_RINGS_H
#define LEMONADE_UPDATER_CREATE_TWO_CONCATENATED_RINGS_H
/**
 * @file
 *
 * @class UpdaterCreateTwoConcatenatedRings
 *
 * @brief create Updater for creation of two concatenated ring polymers
 * 
 * @details A single ring is created in a cubix box of arbitrary size togehter with a second concatenated ring polymer of arbitrary 
 * length 
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>

#include <vector>

template<class IngredientsType>
class UpdaterCreateTwoConcatenatedRings: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
  
public:
  UpdaterCreateTwoConcatenatedRings(IngredientsType& ingredients_, uint32_t numMonomerInRingOne_, uint32_t numMonomerInRingTwo_, uint32_t box_x, uint32_t box_y, uint32_t box_z);
  
  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();
  
  bool checkBondVectorBetweenTwoMonomers(uint32_t indexA, uint32_t indexB);


private:
  using BaseClass::ingredients;
  
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::addMonomerAtPosition;
  using BaseClass::addMonomerInsideConnectedPair;
  using BaseClass::linearizeSystem;
  using BaseClass::moveSystem;

  //! Number of Monomers in first Ring
  uint32_t  numMonomerInRingOne;
  
  //! Number of Monomers in second Ring
  uint32_t  numMonomerInRingTwo;
  
  //! simulation box sizes
  uint32_t boxX,boxY,boxZ;
  
};

/** 
* @brief Constructor handling the new systems paramters
*
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param L_ linear chain length
* @param box_ size of the box
*/
template < class IngredientsType >
UpdaterCreateTwoConcatenatedRings<IngredientsType>::UpdaterCreateTwoConcatenatedRings(IngredientsType& ingredients_, uint32_t numMonomerInRingOne_, uint32_t numMonomerInRingTwo_,  uint32_t box_x, uint32_t box_y, uint32_t box_z):
BaseClass(ingredients_), numMonomerInRingOne(numMonomerInRingOne_), numMonomerInRingTwo(numMonomerInRingTwo_), boxX(box_x), boxY(box_y), boxZ(box_z)
{}

/**
* The initialize function handles the new systems information.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateTwoConcatenatedRings<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterCreateTwoConcatenatedRings" << std::endl;
  
  //set box size
  ingredients.setBoxX(boxX);
  ingredients.setBoxY(boxY);
  ingredients.setBoxZ(boxZ);
  
  //set periodicity
  ingredients.setPeriodicX(false);
  ingredients.setPeriodicY(false);
  ingredients.setPeriodicZ(false);
  
  ingredients.modifyMolecules().setAge(0);
  //add Bondset
  // supress std::cout output of addBondset
  std::streambuf *old = std::cout.rdbuf(); // <-- save
  std::stringstream ss;
  std::cout.rdbuf (ss.rdbuf());       // <-- redirect
  ingredients.modifyBondset().addBFMclassicBondset();
  std::cout.rdbuf (old);
  
  ingredients.synchronize(ingredients);
  
  

 // execute();
}

/**
* Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterCreateTwoConcatenatedRings<IngredientsType>::execute(){
  std::cout << "execute UpdaterCreateTwoConcatenatedRings" << std::endl;
  
  //check if system was already created
  if(ingredients.getMolecules().size()!=0)
    return true;
    
  // check if the number of monomers is at least 8 in every ring
  if(numMonomerInRingOne < 8)
        return true;
 
  if(numMonomerInRingTwo < 8)
        return true;
   
  // create first initial ring of 8 monomers
  addMonomerAtPosition(VectorInt3(boxX/2  , boxY/2  , boxZ/2  ), 0);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2  , boxZ/2  ), 0); ingredients.modifyMolecules().connect( 0, 1);
  addMonomerAtPosition(VectorInt3(boxX/2+4, boxY/2  , boxZ/2  ), 0); ingredients.modifyMolecules().connect( 1, 2);
  addMonomerAtPosition(VectorInt3(boxX/2+4, boxY/2  , boxZ/2-2), 0); ingredients.modifyMolecules().connect( 2, 3);
  
  addMonomerAtPosition(VectorInt3(boxX/2+4, boxY/2  , boxZ/2-4), 0); ingredients.modifyMolecules().connect( 3, 4);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2  , boxZ/2-4), 0); ingredients.modifyMolecules().connect( 4, 5);
  addMonomerAtPosition(VectorInt3(boxX/2  , boxY/2  , boxZ/2-4), 0); ingredients.modifyMolecules().connect( 5, 6);
  addMonomerAtPosition(VectorInt3(boxX/2  , boxY/2  , boxZ/2-2), 0); ingredients.modifyMolecules().connect( 6, 7); ingredients.modifyMolecules().connect( 7, 0);
    
  // create second initial concatenated ring of 8 monomers
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2-2, boxZ/2-2), 0);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2  , boxZ/2-2), 0); ingredients.modifyMolecules().connect( 8, 9);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2+2, boxZ/2-2), 0); ingredients.modifyMolecules().connect( 9, 10);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2+2, boxZ/2-4), 0); ingredients.modifyMolecules().connect( 10, 11);
  
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2+2, boxZ/2-6), 0); ingredients.modifyMolecules().connect( 11, 12);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2  , boxZ/2-6), 0); ingredients.modifyMolecules().connect( 12, 13);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2-2, boxZ/2-6), 0); ingredients.modifyMolecules().connect( 13, 14);
  addMonomerAtPosition(VectorInt3(boxX/2+2, boxY/2-2, boxZ/2-4), 0); ingredients.modifyMolecules().connect( 14, 15); ingredients.modifyMolecules().connect( 15, 8);
  
  
  
  std::vector<int> idxStartMonomerRing;
  idxStartMonomerRing.push_back(0);
  idxStartMonomerRing.push_back(8);
  
  std::vector<int> idxLastMonomerRing;
  idxLastMonomerRing.push_back(7);
  idxLastMonomerRing.push_back(15);


  // successively adding monomers of the first ring between first and last monomer
    for(int numMonos = 0; numMonos < (numMonomerInRingOne-8); numMonos++)
    {
        std::cout << "Adding Mono " << (numMonos+9) << " of " << numMonomerInRingOne << std::endl;
        
        // check for bond and connect
        int idxLast = idxLastMonomerRing.at(0);
        int idxStart = idxStartMonomerRing.at(0);
        
        int tag = (ingredients.getMolecules()[idxLast].getAttributeTag()%2)+1;
        
        while(!addMonomerInsideConnectedPair(idxStart, idxLast,tag))
        {
            //std::cout << "Searching for BondVector between ZZ-DD:" << idxZZ << "-" << idxDD << std::endl;
            moveSystem(1);
        }
                              
        int idx = ingredients.getMolecules().size()-1;
        idxLastMonomerRing.at(0)=idx;
                                                                
    }
      
    // successively adding monomers of the second ring between first and last monomer
    for(int numMonos = 0; numMonos < (numMonomerInRingTwo-8); numMonos++)
    {
        std::cout << "Adding Mono " << (numMonos+9) << " of " << numMonomerInRingTwo << std::endl;
                
        // check for bond and connect
        int idxLast = idxLastMonomerRing.at(1);
        int idxStart = idxStartMonomerRing.at(1);
                
        int tag = (ingredients.getMolecules()[idxLast].getAttributeTag()%2)+1;
                
        while(!addMonomerInsideConnectedPair(idxStart, idxLast,tag))
        {
            //std::cout << "Searching for BondVector between ZZ-DD:" << idxZZ << "-" << idxDD << std::endl;
            moveSystem(1);
        }
                                              
        int idx = ingredients.getMolecules().size()-1;
        idxLastMonomerRing.at(1)=idx;
                                                      
    }
    


  std::cout << "Finial sync...";
  ingredients.synchronize();
  std::cout << "done." << std::endl;

  std::cout << "Finial linearization...";
  linearizeSystem();

  
  for(int i = 0; i < ingredients.getMolecules().size(); i++)
 	  {
	  	  
        ingredients.modifyMolecules()[i].setAttributeTag(0);
	  	  
 	  }
 	  
 
  return true; 
}

/**
* Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateTwoConcatenatedRings<IngredientsType>::cleanup(){
  
}

template < class IngredientsType >
bool UpdaterCreateTwoConcatenatedRings<IngredientsType>::checkBondVectorBetweenTwoMonomers(uint32_t indexA, uint32_t indexB){

	//check the new bondvector bewten the new monomer and indexB
		VectorInt3 checkBV(ingredients.getMolecules()[indexA]-ingredients.getMolecules()[indexB]);
		if( (checkBV.getLength() < 3) && (ingredients.getBondset().isValidStrongCheck(checkBV)) )
			return true;

		return false;
}

#endif /* LEMONADE_UPDATER_CREATE_TWO_CONCATENATED_RINGS_H */
