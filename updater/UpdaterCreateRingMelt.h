/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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

#ifndef LEMONADE_UPDATER_CREATE_RING_MELT_H
#define LEMONADE_UPDATER_CREATE_RING_MELT_H
/**
 * @file
 *
 * @class UpdaterCreateRingMelt
 *
 * @brief create Updater for a single Dendrimer in a melt of linear chains
 * 
 * @details A single dendrimer with arbitrary generation, spacer length and functionality 
 * is created in a cubix box of arbitrary size togehter with linear chains of arbitrary 
 * length with arbitrary concentration
 *
 * @tparam IngredientsType
 *
 **/

#include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>

#include <vector>

template<class IngredientsType>
class UpdaterCreateRingMelt: public UpdaterAbstractCreate<IngredientsType>
{
  typedef UpdaterAbstractCreate<IngredientsType> BaseClass;
  
public:
  UpdaterCreateRingMelt(IngredientsType& ingredients_, uint32_t numRings_, uint32_t numMonomerInRing_, uint32_t box_x, uint32_t box_y, uint32_t box_z);
  
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

  //! Number of Ring Polymers
  uint32_t numRings;

  //! Number of Monomers in Ring
  uint32_t  numMonomerInRing;
  
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
UpdaterCreateRingMelt<IngredientsType>::UpdaterCreateRingMelt(IngredientsType& ingredients_, uint32_t numRings_, uint32_t numMonomerInRing_,  uint32_t box_x, uint32_t box_y, uint32_t box_z):
BaseClass(ingredients_), numRings(numRings_), numMonomerInRing(numMonomerInRing_), boxX(box_x), boxY(box_y), boxZ(box_z)
{}

/**
* The initialize function handles the new systems information.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateRingMelt<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterCreateRingMelt" << std::endl;
  
  //set box size
  ingredients.setBoxX(boxX);
  ingredients.setBoxY(boxY);
  ingredients.setBoxZ(boxZ);
  
  //set periodicity
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  
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
bool UpdaterCreateRingMelt<IngredientsType>::execute(){
  std::cout << "execute UpdaterCreateRingMelt" << std::endl;
  
  //check if system was already created
  if(ingredients.getMolecules().size()!=0)
    return true;
  
  // create first initial rings of three monomers
  std::vector<int> idxStartMonomerRing;
  std::vector<int> idxLastMonomerRing;

  for(int i = 0; i < numRings; i++)
  {
	  int idx = 0;
	  // add start ring
	  while( !addSingleMonomer(1))
	  {
	      moveSystem(1);
	  }
	  idx = ingredients.getMolecules().size()-1;
	  idxStartMonomerRing.push_back(idx);

	  // add additional monomer
	  while( !addMonomerToParent(idx, 2))
	  {
	  	moveSystem(1);
	  }
	  idx = ingredients.getMolecules().size()-1;

	  // add last monomer and close ring
	  while( !addMonomerToParent(idx, 1))
	  {
	 	  	moveSystem(1);
	  }
	  idx = ingredients.getMolecules().size()-1;
	  idxLastMonomerRing.push_back(idx);

	  // check for bond and connect
	  int idxLast = idxLastMonomerRing.at(i);
	  int idxStart = idxStartMonomerRing.at(i);

	  while(!checkBondVectorBetweenTwoMonomers(idxStart,idxLast))
	  {
		  //std::cout << "Searching for BondVector between ZZ-DD:" << idxZZ << "-" << idxDD << std::endl;
		  moveSystem(1);
	  }
	  std::cout << "Ring " << i << " : Found BondVector between Start-Last:" << (ingredients.getMolecules()[idxStart]-ingredients.getMolecules()[idxLast]) << std::endl;

	  ingredients.modifyMolecules().connect( idxStart,idxLast);
  }


  // successively adding monomers between first and last monomer
  for(int numMonos = 0; numMonos < (numMonomerInRing-3); numMonos++)
  {
	  std::cout << "Adding Mono " << (numMonos+4) << " of " << numMonomerInRing << std::endl;

	  for(int i = 0; i < numRings; i++)
	  {
		  	  // check for bond and connect
		  	  int idxLast = idxLastMonomerRing.at(i);
		  	  int idxStart = idxStartMonomerRing.at(i);

		  	  int tag = (ingredients.getMolecules()[idxLast].getAttributeTag()%2)+1;

		  	  while(!addMonomerInsideConnectedPair(idxStart, idxLast,tag))
		  	  {
		  		  //std::cout << "Searching for BondVector between ZZ-DD:" << idxZZ << "-" << idxDD << std::endl;
		  		  moveSystem(1);
		  	  }
		  	  int idx = ingredients.getMolecules().size()-1;
		  	  idxLastMonomerRing.at(i)=idx;
	  }
  }

  std::cout << "Finial sync...";
  ingredients.synchronize();
  std::cout << "done." << std::endl;

  std::cout << "Finial linearization...";
  linearizeSystem();

  for(int i = 0; i < numRings; i++)
 	  {
	  	  for(int numMonos = 0; numMonos < numMonomerInRing; numMonos++)
	  	  {
	  		ingredients.modifyMolecules()[i*numMonomerInRing+numMonos].setAttributeTag(0);
	  	  }
 	  }
/*
  for(int i = 0; i < numRings; i++)
  	  {
	  	  ingredients.modifyMolecules()[i*numMonomerInRing+28-1].setAttributeTag(2);
	  	  ingredients.modifyMolecules()[i*numMonomerInRing+28+26-1].setAttributeTag(2);
	  	  ingredients.modifyMolecules()[i*numMonomerInRing+28+26+26+28-1].setAttributeTag(2);
	  	  ingredients.modifyMolecules()[i*numMonomerInRing+28+26+26+28+26-1].setAttributeTag(2);
  	  }
  */
  //ingredients.synchronize();

  /*
  for(int z = 0; z < ingredients.getMolecules().size(); z++)
  {
	  ingredients.modifyMolecules()[z].setAttributeTag(0);
  }

  for(int z = 0; z < ingredients.getMolecules().size(); z++)
    {
	  if(ingredients.getMolecules().getNumLinks(z) == 1)
  	  ingredients.modifyMolecules()[z].setAttributeTag(2);
    }


  // adding cross-linker to the system
  for(int cl=0; cl < crosslinker; cl++)
  {
	  // Cross-linker are B-Types
	  while( !addSingleMonomer(1))
	  {
		  moveSystem(1);
	  }
	  std::cout << "added cross-linker :" <<  (ingredients.getMolecules().size()-1)  << std::endl;

  }


  std::cout << "Creation done with " << ingredients.getMolecules().size() << " monomers" << std::endl;

  std::cout << "Finial sync...";
  ingredients.synchronize();
  std::cout << "done." << std::endl;

  std::cout << "Finial linearization...";

  linearizeSystem();


  std::cout << "done." << std::endl;
 */
  return true;
}

/**
* Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterCreateRingMelt<IngredientsType>::cleanup(){
  
}

template < class IngredientsType >
bool UpdaterCreateRingMelt<IngredientsType>::checkBondVectorBetweenTwoMonomers(uint32_t indexA, uint32_t indexB){

	//check the new bondvector bewten the new monomer and indexB
		VectorInt3 checkBV(ingredients.getMolecules()[indexA]-ingredients.getMolecules()[indexB]);
		if( (checkBV.getLength() < 3) && (ingredients.getBondset().isValidStrongCheck(checkBV)) )
			return true;

		return false;
}

#endif /* LEMONADE_UPDATER_CREATE_RING_MELT_H */
