#ifndef AnalyzerScatteringSingleObject_H_
#define AnalyzerScatteringSingleObject_H_

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

/*****************************************************************************/
/**
 * @file
 * @brief Analyzer for static struture factor -> Scattering Intensity
 * @author Martin Wengenmayr, Ron Dockhorn
 * @date 24.03.2015, 17.07.2024
 * */
/*****************************************************************************/

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <complex>

/*****************************************************************************/
/**
 * @class AnalyzerScatteringSingleObject
 * @brief  Analyzer class for single objects
 * @details Calculates formfactor via isotrop scattering simulation
 * */
/*****************************************************************************/
template<class IngredientsType>
class AnalyzerScatteringSingleObject: public AbstractAnalyzer {
public:
  //constructor
  AnalyzerScatteringSingleObject(const IngredientsType&, std::string filename = "ScatterFactor.dat", long evalulation_time_=0, int32_t relaxtime_=0, double binWidth_=1);
  
  //destructor
  ~AnalyzerScatteringSingleObject(){}
  
  //open files
  void initialize();
  
  //perform calculations
  bool execute();
  
  //cleanup
  void cleanup();

private:
  const IngredientsType& ingredients;
  /*const typename IngredientsType::molecules_type& molecules;
  const IngredientsType& getIngredients() const {return ingredients;}
  const typename IngredientsType::molecules_type& getMolecules() const {return ingredients.getMolecules();}*/

  // typedef for complex numbers
  typedef std::complex<double> double_complex;

  //! name of output file containing scattering factor results
  std::string outputFile;
  //! default file name if none is assigned during initialization
  // const std::string default_outfilename = "ScatterFactor.dat";
  
  // Random Number Generator (RNG)
  RandomNumberGenerators rng;
  
  int32_t currentTimestep;
  int32_t relaxtime;
  double binWidth;
  uint32_t num_of_q;
  
  // vector for averaged squared absolute value of scattering amplitude
  std::vector<double> averagedSquaredAbsC_q;
  
  // vector for averaged squared absolute value of scattering amplitude
  std::vector<double> averagedSquaredAbsC_q_elements;

  // vector for q-vector absolute values
  std::vector<double> q_factor;
  
  // private functions
  void CalcScatteringAmplitude();
  void CalcBinning();
  void Init_qfactor();
  VectorDouble3 GetRandQVector();
  

  long evalulation_time;
};

/******************************************************************************/
/**
 * @fn void AnalyzerScatteringSingleObject<IngredientsType>::AnalyzerScatteringSingleObject()
 * @brief Constructor: declaration of intern molecule properties
 */
template<class IngredientsType>
AnalyzerScatteringSingleObject<IngredientsType>::AnalyzerScatteringSingleObject(
  const IngredientsType& ingredients_, std::string filename, long evalulation_time_,  int32_t relaxtime_, double binWidth_):
  ingredients(ingredients_), 
  currentTimestep(0),
  relaxtime(relaxtime_),
  binWidth(binWidth_),
  num_of_q(500),
  q_factor(),
  averagedSquaredAbsC_q(num_of_q,0.0),
  averagedSquaredAbsC_q_elements(num_of_q,0.0),
  evalulation_time(evalulation_time_),
  outputFile(filename)
{}
  
/******************************************************************************/
/**
 * @fn void AnalyzerScatteringSingleObject<IngredientsType>::initialize()
 * @brief preparation of output files and memory allocation
 */
template<class IngredientsType>
void AnalyzerScatteringSingleObject<IngredientsType>::initialize(){
  std::cout << "\nAnalyzerScatteringSingleObject initialise\n";
  
  Init_qfactor();

  std::cout << "AnalyzerScatteringSingleObject initialised successfully\n\n";
  //first config is written in by BFM file reader in initialise, so execute has to be called in this step the first time too
  execute();
}

/******************************************************************************/
/**
 * @fn void AnalyzerScatteringSingleObject<IngredientsType>::execute()
 * @brief Calculates and writes out the Rg^2 and Ree^2 
 */
template<class IngredientsType>
bool AnalyzerScatteringSingleObject<IngredientsType>::execute(){
  currentTimestep++;
  std::cout << "AnalyzerScatteringSingleObject starts to execute timestep nr "<<currentTimestep <<"\n";
  
  if(ingredients.getMolecules().getAge() > evalulation_time)
  {
    CalcScatteringAmplitude();
  }
  

  std::cout << "AnalyzerScatteringSingleObject executed successfully timestep nr "<<currentTimestep
  << " with age " << ingredients.getMolecules().getAge()<<std::endl<<std::endl;
  
  return true;
}

/******************************************************************************/
/**
 * @fn void AnalyzerScatteringSingleObject<IngredientsType>::cleanup()
 * @brief Write out results
 */
template<class IngredientsType>
void AnalyzerScatteringSingleObject<IngredientsType>::cleanup(){
  
  std::cout << "AnalyzerScatteringSingleObject starts clean up after timestep "<<currentTimestep <<"\n";
  std::ofstream FormFactorFile;
  
  if(currentTimestep > relaxtime){

	// get the filename and path
	// std::string filenameGeneral=ingredients.getName();
	// // delete the .bfm in the name
	// filenameGeneral.erase (ingredients.getName().length()-4, ingredients.getName().length());
 //
	// //new filename
	// std::string filename_ScatteringFct = filenameGeneral + "_ScatteringFct.dat";

  u_int32_t numScatteringObj=0;

  for(uint32_t k=0;k<ingredients.getMolecules().size();k++){
    //if((ingredients.getMolecules()[k].getAttributeTag()==1))
     numScatteringObj++;
  }

    FormFactorFile.open (outputFile.c_str(),std::ios::out);
    FormFactorFile << "# Molecular Scattering Function\n"
    << "# q    S(q)   samples"<< std::endl;
    for (uint32_t i=0;i<q_factor.size();i++){
      FormFactorFile << q_factor.at(i)*((2*M_PI)/ingredients.getBoxX()) <<" "
      << averagedSquaredAbsC_q.at(i)/(averagedSquaredAbsC_q_elements.at(i)*(1.0*numScatteringObj)) << " "
      << averagedSquaredAbsC_q_elements.at(i)
	  << std::endl;
    }
  }else{
    std::cout<<"****   **** No output file written!! ****   ****\n****   **** Not enought configs in input file!! ****   ****\n";
  }

  std::cout << "AnalyzerScatteringSingleObject cleanup successfully performed\n\n";
}

/******************************************************************************/
/**
 * @fn void AnalyzerScatteringSingleObject<IngredientsType>::Init_qfactor( void )
 * @brief initialise set of q vectors in multiplicative manner
 */
template<class IngredientsType>
void AnalyzerScatteringSingleObject<IngredientsType>::Init_qfactor(){
  double k(std::exp(std::log((10.0*ingredients.getBoxX()))*1.0/(num_of_q)));
  double x(0.1);
  std::cout << "multiplicator for abs(q) k = " << k << std::endl;
  for(uint32_t i=0; i<num_of_q; i++){
    x*=k;
    q_factor.push_back(x);
    std::cout << "multiplicator for abs(q) x = " << x << std::endl;
  }
  if((q_factor.at(num_of_q-1) - (double)(4*ingredients.getBoxX())) > 0.00001)
    throw std::runtime_error("wrong q_vector array!");
}

/******************************************************************************/
/**
 * @fn void AnalyzerScatteringSingleObject<IngredientsType>::CalcScatteringAmplitude( void )
 * @brief adding a count per strand length interval 
 */
template<class IngredientsType>
void AnalyzerScatteringSingleObject<IngredientsType>::CalcScatteringAmplitude(){
	// loop over some randomly catched q vectors
	// for(uint32_t i=0;i<num_of_q;i++){
	for(uint32_t i=0;i<50;i++){

		// get random unit q-vector
		VectorDouble3 q(GetRandQVector());

		//loop over some different absolute values of q
		for(uint32_t j=0;j<q_factor.size();j++){
			// multiply q with a certain number to get equidistant points in log-log plot
			VectorDouble3 q_j=(((2.0*M_PI)/ingredients.getBoxX())*q_factor.at(j))*q;
			// scattering amplitude of C_q
			double_complex C_q(0.0,0.0);
			//loop over molecules
			for(uint32_t k=0;k<ingredients.getMolecules().size();k++){
				// "cast" position vector to double vector
				VectorDouble3 r (ingredients.getMolecules()[k]);
				//scalar product of q and monomer position
				double phi( (double)((-1)*(r*q_j)) );
				//"cast" to imaginary complex number
				double_complex iphi(0.0,phi);

				// get factor:
				double factor = 1.0;
				//if((ingredients.getMolecules()[k].getAttributeTag()==1))
					//factor = 1.0;

				//if((ingredients.getMolecules()[k].getAttributeTag()!=1))
					//factor = 0.0;

				// sum over the squared abs of the exponential of iphi
				C_q+=factor*std::exp(iphi);
			} /* end loop over molecule */
			// add C_q depending on q as part of the scattering function to the scattering function container
			averagedSquaredAbsC_q.at(j)+=((double)(std::real(C_q*std::conj(C_q))));
			averagedSquaredAbsC_q_elements.at(j)+=1.0;
		}/* end loop over differnt absolute values of q */

	}/* end loop over random q vectors */
}

/******************************************************************************/
/**
 * @fn VectorDouble3 AnalyzerScatteringSingleObject<IngredientsType>::GetRandQVector()
 * @brief get a randomly orientated 3D Vector 
 * @return VectorDouble3 q, restricted to the unit sphere
 */
template<class IngredientsType>
VectorDouble3 AnalyzerScatteringSingleObject<IngredientsType>::GetRandQVector(){
  VectorDouble3 q(1.0,1.0,1.0);
  //get three random numbers in [-1,1] 
  double a((rng.r250_drand()*2)-1);
  double b((rng.r250_drand()*2)-1);
  double c((rng.r250_drand()*2)-1);
  //check them to be inside the unit sphere
  if((a*a+b*b+c*c)<=1){
    // if yes, set them as coordinates for q
    q.setX(a); q.setY(b); q.setZ(c);
  }else if((a*a+b*b+c*c)>1){
    // if not, do an iterative call of GetRandQVector to get another set of random coordinates
    q=GetRandQVector();
  }
  q.normalize();
  return q;
}

#endif /* AnalyzerScatteringSingleObject_H_ */
